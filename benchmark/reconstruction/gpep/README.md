# Focal-recovery two-view geometry benchmarks

End-to-end ETH3D benchmarks for the two dispatch branches of
`EstimateTwoViewGeometry` that recover focal length. Each runs the
reconstruction benchmark twice, once per colmap binary, and prints a
setup check, a runtime comparison and the paired accuracy tables.

| Script | Estimator | Reached when | Scenes |
| --- | --- | --- | --- |
| `bench_shared_focal.py` | `EstimateSharedFocalTwoViewGeometry` | both images share one camera with no focal prior | 19 |
| `bench_one_sided_focal.py` | `EstimateOneSidedFocalTwoViewGeometry` | exactly one of the two cameras has a focal prior | 25 |

`--help` on either lists every option.

## Setup

**No stock benchmark configuration reaches either estimator.** The harness
primes either every camera with its GT focal or none, and the extractor gives
ETH3D one camera per image, so a plain run yields no shared camera and no
mixed pair. Each script therefore wraps the colmap binary in a generated
shell shim (written to `<run-dir>/.shims/`) that adjusts the database between
the harness's extraction and matching phases:

* `bench_shared_focal.py` passes `--single_camera 1` and clears the focal
  prior before matching.
* `bench_one_sided_focal.py` gives a deterministic `--calibrated-fraction`
  (default 0.5) of each scene's cameras their GT intrinsics and resets the
  rest to the extractor's default focal.

Nothing under `benchmark/` is modified. The shim also times each phase into
`<run-dir>/timing-<arm>.csv`.

**Scene selection.** One camera per image is a property of the extractor, not
of a scene, so `bench_one_sided_focal.py` needs no selection and runs all 25
DSLR scenes; its split is per image, so scenes with several GT cameras are
fine. `bench_shared_focal.py` forces `--single_camera`, which is only true of
the 19 scenes whose GT has exactly one camera, so it restricts to those:

```
botanical_garden  boulders    bridge      delivery_area  door
exhibition_hall   kicker      lecture_room  living_room   lounge
meadow            observatory office       old_computer   pipes
relief            relief_2    statue       terrace_2
```

The six excluded have 2 to 6 GT cameras: `courtyard`, `electro`, `facade`,
`playground`, `terrace`, `terrains`. The script derives this list from the GT
at startup rather than hard-coding it, so it cannot go stale; pass `--scenes`
to override.

### Prerequisites

Commands below are run from this directory.

```bash
python ../download.py --datasets eth3d   # DSLR scenes are ~11 GB extracted
pip install -r ../requirements.txt       # needs pycolmap
```

Two colmap binaries: **with** the solvers (stock main) and **without**.

Both solvers are already merged, so the without-binary is not an older commit
-- the last commit predating both is far enough back that unrelated work
would land in the measured difference. Build it from current main with the
two dispatch branches disabled. In
`src/colmap/estimators/two_view_geometry.cc`, add to the anonymous namespace:

```c++
constexpr bool kEnableFocalSolvers = false;
```

and gate both branches on it:

```c++
    if (kEnableFocalSolvers &&
        IsCameraCalibrated(camera1) != IsCameraCalibrated(camera2) && ...

    } else if (kEnableFocalSolvers &&
               camera1.camera_id == camera2.camera_id && ...
```

Both configurations then fall through to
`EstimateUncalibratedTwoViewGeometry`, the plain fundamental-matrix path.
This affects pinhole pairs only: a mixed spherical/pinhole pair still reaches
`EstimateSphericalTwoViewGeometry`.

That build is `--baseline-colmap`, your own is `--branch-colmap`. Use the
same cmake options for both.

### Running

`--mapper` selects the pipeline; it defaults to `incremental`. Give each
mapper its own `--run-name`, since they cannot share a workspace.

```bash
# Smoke test first: one scene, one seed
python bench_shared_focal.py --baseline-colmap A --branch-colmap B \
    --scenes door --num-seeds 1 --run-name smoke

# Incremental, all scenes, five seeds
python bench_shared_focal.py --baseline-colmap A --branch-colmap B \
    --mapper incremental --num-seeds 5 --run-name shared-incremental

# Global
python bench_shared_focal.py --baseline-colmap A --branch-colmap B \
    --mapper global --num-seeds 5 --run-name shared-global
```

`bench_one_sided_focal.py` takes the same arguments, plus
`--calibrated-fraction` (default 0.5) for the share of cameras given a GT
focal.

Both arms share `--run-path/--run-name`, so features and raw matches are
extracted once, and they run sequentially. Each arm runs seeds `0..N-1` and
the difference is taken per seed before averaging. The results below are all
`--mapper incremental`.

## Results

Both mappers, measured on one RTX 3080 at `--num-parallel-scenes 1`, five seeds, both arms
from one worktree at `af20f0fa` differing only in `kEnableFocalSolvers`.
`nofocal` is the without-arm, `focal` is stock main, and `A - B` is
`focal - nofocal`, so positive favours the solver.

### Runtime

| mapper | | verification, without | verification, with | change |
| --- | --- | --- | --- | --- |
| incremental | shared focal | 9.9 ms/pair | 95.7 ms/pair | +862% |
| incremental | one-sided focal | 18.7 ms/pair | 46.8 ms/pair | +150% |
| global | shared focal | 16.5 ms/pair | 97.9 ms/pair | +494% |
| global | one-sided focal | 29.9 ms/pair | 56.2 ms/pair | +88% |

The denominator is the pairs verification ran on, those holding raw matches:
10662 and 15422 incremental, 7848 and 11670 global. The incremental figures
are one seed reusing the accuracy run's workspaces; the global ones are the
mean over all five seeds of their own runs, which varied by under 2%.

The incremental figure is a matching-phase total that at `--quality high`
also includes guided matching; the global pipeline disables guided matching,
so its figure is verification alone.

The two mappers' absolute values are not comparable, for two reasons that
compound.

Verification itself is slower under the global option set. Re-running one
scene (`exhibition_hall`) from the same workspace with the same binary, so
the pair set is identical and only the mapper's options differ, takes 35.6 s
with the incremental options and 48.8 s with the global ones -- and the
global run skips guided matching entirely, so the estimation-only gap is
wider than that 37%. The likely driver is `max_error`, 4.0 px against 1.0 px:
a tighter threshold lowers the inlier ratio, and RANSAC then needs more
trials to reach the same confidence.

The denominator also differs. `min_num_inliers` is 15 for incremental and 30
for global, and a pair below it keeps no matches: the incremental workspaces
hold 10662 pairs with matches of which 2814 have fewer than 30, and the
global ones hold exactly the 7848 that remain. The dropped pairs are the
weakest and cheapest, so removing them raises the mean.

Isolating one scene (`kicker`, 465 pairs) gives 5.6 s without against 63.3 s
with, for 276 against 263 verified pairs and 214k against 216k inliers, so
the difference is in estimation rather than in downstream matching.

The two baselines are different workloads, so compare within a row rather
than down a column: the shared-focal setting is 100% uncalibrated pairs,
while the one-sided setting is 24.9% both-calibrated (essential matrix in
both arms), 51.4% mixed and 23.7% uncalibrated, over a different scene set.

### bench_shared_focal.py

Setup check -- one camera and no prior per scene:

```
scene                 images  cameras  priors    pairs
botanical_garden          30        1       0      137
boulders                  26        1       0      178
bridge                   110        1       0     1489
...
(all 19 scenes: 1 camera, 0 priors, non-zero pairs)
```

**Incremental mapper:**

```
B = nofocal  (AUC mean ± std over 5 seeds)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  20.01± 0.59  31.40± 0.81  46.91± 1.37  49.21± 1.49
boulders           5  25.97± 0.15  50.79± 0.20  85.66± 0.07  91.21± 0.04
bridge             5  53.18± 0.16  75.49± 0.09  94.92± 0.03  97.46± 0.02
delivery_area      5  33.82± 0.15  56.85± 0.12  89.46± 0.11  94.33± 0.15
door               5  28.84± 4.83  57.62± 3.12  87.63± 0.74  92.52± 1.35
exhibition_hall    5   7.94± 0.47  25.21± 0.77  77.65± 0.31  87.71± 0.20
kicker             5  28.75± 2.25  50.25± 3.66  81.68± 1.39  87.32± 0.72
lecture_room       5  41.34± 3.75  60.06± 3.64  87.53± 2.29  93.01± 1.74
living_room        5  44.95± 1.93  68.77± 1.13  92.73± 0.28  96.31± 0.15
lounge             5  29.62± 2.43  37.46± 2.36  43.76± 2.82  44.55± 2.90
meadow             5  21.16± 5.61  41.63± 8.95  78.23± 9.31  86.30± 9.61
observatory        5  21.05± 0.20  41.21± 0.29  83.87± 0.14  91.74± 0.08
office             5  12.39± 0.34  25.53± 0.61  69.94± 1.29  83.09± 1.24
old_computer       5  17.19± 0.57  35.12± 0.59  80.32± 0.62  89.54± 0.48
pipes              5  31.45± 0.57  57.70± 0.49  89.99± 0.17  95.03± 0.08
relief             5  34.47± 0.36  62.43± 0.29  91.62± 0.07  95.81± 0.04
relief_2           5  30.95± 0.65  60.49± 0.53  91.70± 0.13  95.85± 0.06
statue             5  33.78± 0.20  62.97± 0.21  92.59± 0.04  96.30± 0.02
terrace_2          5  70.17± 0.35  83.01± 0.22  96.27± 0.05  98.14± 0.02
------------------------------------------------------------------------
average            5  30.90± 0.33  51.79± 0.45  82.24± 0.39  87.65± 0.37
```

```
A - B = focal - nofocal  (AUC mean ± std over 5 seeds, shared seeds -> paired)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  +0.31± 1.28  +0.23± 1.82  -0.16± 2.53  -0.23± 2.64
boulders           5  +0.62± 0.26  +1.28± 0.30  +2.61± 0.09  +2.79± 0.07
bridge             5  -0.07± 0.45  +0.01± 0.21  +0.05± 0.05  +0.03± 0.03
delivery_area      5  -0.25± 0.43  -0.29± 0.75  -0.50± 1.07  -0.55± 1.11
door               5  -3.74± 5.74  -2.52± 3.46  -0.60± 0.85  -1.39± 1.35
exhibition_hall    5  -0.09± 0.63  -0.23± 1.25  -0.01± 0.55  +0.04± 0.29
kicker             5  -0.34± 2.76  -0.69± 3.74  -0.28± 1.50  -0.24± 0.91
lecture_room       5  -0.59± 3.95  -0.20± 3.94  +0.53± 2.47  +0.54± 1.81
living_room        5  -0.26± 1.46  -0.05± 0.97  +0.03± 0.27  +0.01± 0.13
lounge             5  +0.33± 2.53  -0.42± 1.80  -0.21± 0.59  -0.11± 0.30
meadow             5  +7.97± 5.77 +11.13±10.12  +9.72±10.24  +7.63±10.13
observatory        5  -0.17± 0.25  -0.32± 0.35  -0.14± 0.16  -0.07± 0.09
office             5  -0.10± 0.64  -0.37± 0.94  -0.20± 1.24  -0.18± 1.14
old_computer       5  -0.60± 1.08  +0.07± 1.59  +0.46± 1.03  +0.32± 0.65
pipes              5  +1.69± 0.80  +1.15± 0.73  +0.37± 0.28  +0.17± 0.12
relief             5  +0.14± 0.36  +0.07± 0.27  +0.03± 0.07  +0.01± 0.03
relief_2           5  -1.05± 0.83  -0.84± 0.67  -0.20± 0.18  -0.10± 0.09
statue             5  +0.05± 0.14  -0.00± 0.15  -0.00± 0.03  -0.00± 0.01
terrace_2          5  +0.65± 1.16  +0.46± 0.73  +0.11± 0.15  +0.05± 0.08
------------------------------------------------------------------------
average            5  +0.24± 0.44  +0.45± 0.71  +0.61± 0.58  +0.46± 0.51
```

**Global mapper**, same 19 scenes:

```
B = nofocal  (AUC mean ± std over 5 seeds)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  22.75± 3.22  44.50± 8.00  84.32± 4.42  91.47± 2.79
boulders           5  24.71± 0.06  49.32± 0.15  87.73± 0.03  93.67± 0.02
bridge             5  52.09± 0.03  74.83± 0.02  94.76± 0.01  97.37± 0.00
delivery_area      5  27.44± 9.17  46.68±16.02  76.83±27.49  81.60±28.89
door               5  61.11± 0.53  76.83± 0.34  91.56± 0.07  96.16± 0.07
exhibition_hall    5   7.92± 0.08  26.80± 0.09  78.60± 0.02  88.23± 0.01
kicker             5  30.04± 0.22  51.83± 0.21  82.31± 0.05  87.65± 0.03
lecture_room       5  36.89± 0.79  56.31± 0.61  86.70± 0.13  92.82± 0.07
living_room        5  42.12± 0.60  66.65± 0.42  92.11± 0.12  95.96± 0.06
lounge             5  27.19± 1.31  36.82± 0.70  44.70± 0.14  45.68± 0.07
meadow             5  20.75±11.05  40.07±19.88  69.48±29.24  74.80±30.55
observatory        5  22.43± 0.12  42.33± 0.10  84.09± 0.02  91.85± 0.01
office             5   9.86± 0.32  23.01± 0.43  67.85± 0.66  81.94± 0.42
old_computer       5  15.17± 0.16  32.28± 0.21  78.65± 0.14  88.73± 0.08
pipes              5  30.05± 0.22  55.79± 0.23  89.49± 0.07  94.74± 0.05
relief             5  32.54± 0.07  60.10± 0.05  90.93± 0.01  95.46± 0.00
relief_2           5  29.41± 0.07  59.24± 0.04  91.48± 0.00  95.74± 0.00
statue             5  33.80± 0.01  62.95± 0.01  92.59± 0.00  96.30± 0.00
terrace_2          5  71.41± 0.06  83.58± 0.04  96.37± 0.01  98.18± 0.00
------------------------------------------------------------------------
average            5  31.46± 1.00  52.10± 1.97  83.19± 3.02  88.86± 3.14
```

```
A - B = focal - nofocal  (AUC mean ± std over 5 seeds, shared seeds -> paired)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  +0.39± 0.97  -0.74± 3.09  -1.56± 3.14  -1.06± 2.05
boulders           5  +0.01± 0.12  +0.07± 0.12  +0.02± 0.02  +0.00± 0.02
bridge             5  -0.03± 0.04  -0.02± 0.02  -0.00± 0.01  -0.00± 0.00
delivery_area      5  -0.57±15.56  -1.45±27.98  -2.61±48.21  -2.79±50.75
door               5  +0.18± 0.23  +0.07± 0.12  +0.22± 0.48  +0.02± 0.07
exhibition_hall    5  +0.01± 0.08  +0.02± 0.09  +0.00± 0.02  +0.00± 0.01
kicker             5  +0.07± 0.13  +0.04± 0.13  +0.01± 0.03  -0.00± 0.03
lecture_room       5  +0.11± 0.40  +0.01± 0.38  -0.02± 0.15  -0.01± 0.08
living_room        5  -0.59± 1.37  -0.37± 0.91  -0.11± 0.25  -0.06± 0.13
lounge             5  -0.56± 1.52  -0.32± 0.80  -0.06± 0.16  -0.03± 0.08
meadow             5  +0.44±11.52  -2.92±19.85  -4.97±31.80  -3.94±34.22
observatory        5  -0.00± 0.02  +0.05± 0.06  +0.05± 0.06  +0.03± 0.03
office             5  -0.48± 1.06  -1.14± 2.40  -4.89±10.74  -6.04±13.39
old_computer       5  -0.03± 0.14  -0.00± 0.15  +0.01± 0.11  +0.00± 0.06
pipes              5  +0.08± 0.64  +0.08± 0.66  +0.02± 0.16  +0.01± 0.10
relief             5  +0.09± 0.03  +0.11± 0.06  +0.03± 0.02  +0.02± 0.01
relief_2           5  +0.03± 0.11  +0.00± 0.05  -0.00± 0.01  -0.00± 0.00
statue             5  -0.01± 0.01  -0.01± 0.01  -0.00± 0.00  -0.00± 0.00
terrace_2          5  -0.07± 0.12  -0.04± 0.07  -0.01± 0.02  -0.00± 0.01
------------------------------------------------------------------------
average            5  -0.05± 1.42  -0.35± 2.70  -0.73± 4.80  -0.73± 5.16
```

### bench_one_sided_focal.py

Setup check -- `mixed` counts pairs with exactly one calibrated side, the
only pairs the estimator sees:

```
scene                 cameras   calib    pairs    mixed  mixed %
botanical_garden           30      15      137       70    51.1%
boulders                   26      13      178       89    50.0%
bridge                    110      55     1503      750    49.9%
...
TOTAL                                     9324     4752    51.0%
```

**Incremental mapper:**

```
B = nofocal  (AUC mean ± std over 5 seeds)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  29.93± 1.18  39.45± 1.59  48.16± 2.01  49.29± 2.07
boulders           5  53.53± 0.47  75.39± 0.27  94.44± 0.08  97.19± 0.05
bridge             5  60.39± 0.31  79.44± 0.17  95.75± 0.03  97.86± 0.01
courtyard          5  38.60± 1.72  60.35± 1.63  90.69± 0.48  95.32± 0.26
delivery_area      5  43.29± 0.98  63.95± 1.02  91.71± 0.51  95.84± 0.27
door               5  16.46± 1.83  36.54± 1.74  77.09± 0.80  86.17± 0.40
electro            5  35.78± 0.75  59.12± 0.57  89.48± 0.26  94.46± 0.17
exhibition_hall    5  25.90± 1.40  50.19± 1.42  86.51± 0.49  92.46± 0.31
facade             5  50.48± 0.40  70.30± 0.46  93.41± 0.20  96.70± 0.11
kicker             5  38.66± 4.34  59.93± 5.37  84.73± 3.27  89.01± 1.81
lecture_room       5  54.25± 1.45  73.60± 1.13  94.15± 0.50  97.07± 0.25
living_room        5  54.21± 1.89  74.48± 1.75  94.26± 0.68  97.10± 0.38
lounge             5  11.76±11.31  15.61±15.21  19.47±19.47  20.45±19.86
meadow             5   8.57± 6.00  16.98± 9.87  40.10±22.87  48.88±27.19
observatory        5  17.21± 0.31  39.88± 0.84  82.76± 2.12  90.67± 2.32
office             5  19.36± 0.92  32.58± 2.10  61.14± 7.58  70.33± 8.81
old_computer       5  27.12± 1.09  51.35± 1.35  87.91± 0.48  93.81± 0.26
pipes              5  30.56± 4.55  52.95± 7.95  79.40±12.29  83.54±12.43
playground         5  51.63±13.42  69.04±17.31  85.30±20.03  87.68±20.11
relief             5  37.47± 0.57  61.37± 0.62  90.32± 0.18  95.08± 0.09
relief_2           5  42.60± 0.38  65.39± 0.30  91.91± 0.07  95.96± 0.03
statue             5  37.95± 0.32  60.33± 0.25  91.02± 0.09  95.51± 0.05
terrace            5  54.72± 0.83  75.70± 0.49  95.02± 0.07  97.51± 0.04
terrace_2          5  69.53± 1.56  82.34± 1.69  95.76± 1.66  97.51± 1.66
terrains           5  38.00± 3.30  63.74± 2.79  92.13± 0.61  95.98± 0.31
------------------------------------------------------------------------
average            5  37.92± 0.44  57.20± 0.52  82.11± 0.72  86.46± 0.88
```

```
A - B = focal - nofocal  (AUC mean ± std over 5 seeds, shared seeds -> paired)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  +1.48± 1.06  +2.11± 1.35  +2.54± 1.96  +2.62± 2.05
boulders           5  +0.08± 0.27  +0.09± 0.14  +0.07± 0.06  +0.04± 0.04
bridge             5  -1.31± 1.32  -0.94± 0.66  -0.31± 0.24  -0.21± 0.23
courtyard          5  -5.22±14.11  -8.87±21.46 -13.86±31.00 -14.04±31.36
delivery_area      5  -0.50± 0.88  -0.70± 0.94  -0.39± 0.47  -0.20± 0.25
door               5  +0.85± 2.49  -0.27± 1.56  -1.67± 1.08  -0.82± 0.53
electro            5  +0.68± 2.10  +0.82± 2.04  +0.16± 0.59  +0.11± 0.32
exhibition_hall    5  -1.24± 2.17  -0.80± 1.84  -0.26± 0.68  -0.10± 0.42
facade             5  -0.23± 0.47  -0.31± 0.43  -0.08± 0.15  -0.04± 0.08
kicker             5  -0.61± 5.52  +1.00± 5.72  +1.14± 2.86  +0.56± 1.56
lecture_room       5  -0.40± 1.40  -0.30± 1.03  -0.04± 0.45  -0.02± 0.23
living_room        5  +0.73± 3.53  +0.42± 3.16  +0.07± 1.22  +0.03± 0.68
lounge             5  +9.18±10.99 +11.30±14.93 +12.58±19.42 +12.51±19.80
meadow             5  -2.67± 9.42  -1.69±18.21  -2.18±36.61  -2.74±42.67
observatory        5  +0.12± 0.52  +0.12± 0.95  +0.52± 2.13  +0.81± 2.32
office             5  -0.17± 1.99  +0.01± 4.85  -0.31±14.53  -0.34±16.02
old_computer       5  -1.80± 0.86  -2.46± 1.83  -0.96± 0.64  -0.52± 0.32
pipes              5  -4.19± 5.24  -2.75±10.74  +1.21±18.75  +1.59±19.81
playground         5  +4.91±13.53  +7.29±17.29  +9.70±20.03  +9.80±20.12
relief             5  +0.15± 0.47  +0.11± 0.58  +0.04± 0.18  +0.02± 0.09
relief_2           5  +0.18± 0.69  +0.08± 0.59  +0.02± 0.17  +0.01± 0.09
statue             5  -0.80± 1.12  -0.83± 0.65  -0.27± 0.17  -0.14± 0.08
terrace            5  -0.40± 0.96  -0.56± 0.53  -0.18± 0.10  -0.09± 0.05
terrace_2          5  +0.48± 1.58  +0.71± 1.68  +0.75± 1.66  +0.75± 1.66
terrains           5  -0.16± 3.23  -0.54± 3.15  -0.18± 0.83  -0.10± 0.41
------------------------------------------------------------------------
average            5  -0.03± 0.98  +0.12± 1.18  +0.32± 1.19  +0.38± 1.30
```

**Global mapper**, same 25 scenes:


```
B = nofocal  (AUC mean ± std over 5 seeds)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  37.77± 9.20  59.08±15.62  83.40±16.72  87.92±14.72
boulders           5  43.51± 4.35  66.63± 5.63  89.05± 6.62  92.24± 6.68
bridge             5  59.02± 0.21  78.37± 0.19  95.40± 0.10  97.68± 0.08
courtyard          5  33.02± 2.35  55.26± 3.62  87.30± 2.82  92.51± 2.80
delivery_area      5  42.51± 0.21  63.43± 0.24  91.32± 0.07  95.62± 0.04
door               5  25.93± 0.50  44.52± 0.48  79.86± 0.23  87.55± 0.11
electro            5  31.40± 1.28  53.55± 1.57  83.67± 1.70  89.22± 1.91
exhibition_hall    5  25.21± 0.13  49.31± 0.35  85.83± 1.10  91.75± 1.21
facade             5  49.99± 0.19  69.70± 0.57  92.38± 1.95  95.67± 2.14
kicker             5  33.35± 0.59  57.61± 0.78  84.52± 0.18  89.19± 0.19
lecture_room       5  44.83± 6.33  62.35± 8.37  85.74± 8.59  90.02± 8.23
living_room        5  55.04± 0.42  75.19± 0.31  94.52± 0.09  97.26± 0.05
lounge             5  27.79± 3.68  35.59± 4.68  43.22± 3.67  44.81± 2.14
meadow             5   4.20± 2.01   8.91± 5.66  23.50±11.51  33.06±13.45
observatory        5  10.70± 0.12  29.94± 0.38  79.25± 0.19  89.42± 0.10
office             5  17.33± 2.38  33.16± 3.19  74.35± 1.38  85.97± 1.10
old_computer       5  24.19± 0.66  46.93± 0.98  86.51± 0.36  93.08± 0.20
pipes              5  34.22± 3.41  60.95± 5.35  87.61± 7.31  90.95± 7.57
playground         5  50.08± 3.25  71.20± 3.59  91.46± 4.36  94.17± 4.50
relief             5  35.37± 1.23  57.95± 2.27  87.21± 3.27  92.23± 3.41
relief_2           5  36.95± 1.80  58.64± 3.98  85.18± 6.79  89.41± 7.21
statue             5  39.28± 0.07  59.59± 0.13  90.37± 0.08  95.18± 0.04
terrace            5  52.06± 1.64  73.34± 2.59  93.22± 3.62  95.74± 3.76
terrace_2          5  69.52± 0.09  81.98± 0.04  96.16± 0.01  98.08± 0.01
terrains           5  41.10± 1.31  65.02± 1.99  90.13± 3.68  93.58± 3.95
------------------------------------------------------------------------
average            5  36.98± 0.81  56.73± 1.31  83.25± 1.58  88.09± 1.49
```

```
A - B = focal - nofocal  (AUC mean ± std over 5 seeds, shared seeds -> paired)
scene              N     @0.5         @1.0         @5.0        @10.0    
------------------------------------------------------------------------
botanical_garden   5  +4.14± 5.23  +4.43±15.32  +1.92±23.94  +2.35±20.10
boulders           5  +3.69± 4.67  +3.60± 5.84  +3.28± 6.69  +3.16± 6.71
bridge             5  -0.06± 0.38  -0.06± 0.36  -0.04± 0.19  -0.03± 0.15
courtyard          5  +1.74± 2.31  +1.57± 2.97  +2.06± 2.69  +2.11± 2.77
delivery_area      5  +0.11± 0.20  +0.05± 0.22  +0.02± 0.08  +0.01± 0.04
door               5  -0.09± 0.28  -0.02± 0.09  -0.00± 0.02  -0.00± 0.01
electro            5  -0.64± 1.53  -0.48± 2.43  +0.69± 3.79  +1.82± 3.27
exhibition_hall    5  +0.08± 0.14  +0.16± 0.31  +0.50± 1.10  +0.55± 1.20
facade             5  +0.10± 0.14  +0.31± 0.59  +0.89± 1.96  +0.97± 2.15
kicker             5  +0.54± 1.19  +0.48± 1.02  +0.09± 0.22  -0.35± 0.18
lecture_room       5  +1.18± 2.74  +1.32± 3.11  +1.67± 3.66  +1.73± 3.78
living_room        5  +0.11± 0.30  +0.08± 0.22  +0.02± 0.07  +0.01± 0.04
lounge             5  -0.12± 0.66  +0.04± 0.33  +0.01± 0.07  +0.00± 0.03
meadow             5  -1.09± 3.69  -1.47±10.34  -0.49±22.57  -1.93±25.99
observatory        5  -0.01± 0.11  -0.00± 0.25  +0.00± 0.13  +0.00± 0.07
office             5  -1.45± 3.20  -2.15± 4.31  -2.57± 2.72  -2.43± 2.95
old_computer       5  +0.26± 0.87  +0.19± 1.34  +0.09± 0.48  +0.04± 0.26
pipes              5  +0.74± 3.09  +1.77± 4.55  +2.64± 6.01  +2.75± 6.20
playground         5  +0.15± 6.94  -0.17± 8.40  -0.78± 9.57  -0.61± 9.24
relief             5  +0.85± 1.22  +1.63± 2.31  +2.37± 3.29  +2.48± 3.41
relief_2           5  +1.25± 1.53  +2.87± 3.14  +4.60± 4.73  +4.84± 4.96
statue             5  +0.02± 0.05  -0.02± 0.08  -0.01± 0.03  -0.01± 0.02
terrace            5  +0.27± 2.12  +0.15± 4.00  +0.03± 5.71  +0.01± 5.93
terrace_2          5  -0.04± 0.12  -0.07± 0.11  -0.03± 0.04  -0.01± 0.02
terrains           5  +0.71± 1.15  +1.34± 1.95  +2.50± 3.69  +2.67± 3.95
------------------------------------------------------------------------
average            5  +0.50± 0.45  +0.62± 1.16  +0.78± 1.95  +0.80± 1.86
```

## Notes

* `--threads-per-scene 1` is the default and should stay there: RANSAC seeds
  per thread as `random_seed + omp_get_thread_num()`, so only single-threaded
  scenes are reproducible at a fixed seed, which the paired difference relies
  on. Use `--num-parallel-scenes` for throughput.
* `--num-parallel-scenes` defaults to 1. Each concurrent scene holds about
  3.4 GB of GPU memory during matching (two SIFT GPU matcher workers), so
  three scenes cannot fit a 10 GB card. On failure colmap logs `Not enough
  GPU memory to match N features`, the worker never starts, and the scene
  ends with zero verified pairs and registers nothing -- raw matches survive,
  but the two-view geometries just cleared by `--overwrite` are never
  regenerated. The scripts abort rather than print a comparison when any
  scene has zero verified pairs.
* The two mappers measure different things. The incremental mapper consumes
  the focal these solvers estimate (`incremental_mapper_impl.cc`,
  `info.camera1 = two_view_geometry.camera1`), whereas view-graph calibration
  recomputes focals from the view graph and then clears
  `tvg.camera1/camera2` so consumers use its own K. Under the global mapper
  the solvers therefore act only indirectly, through the F they produce and
  the inliers they keep. Global also disables guided matching and uses
  different two-view thresholds (`max_error` 1.0, `min_num_inliers` 30,
  `min_inlier_ratio` 0.25). The raised `min_num_inliers` is why its pair
  counts and baseline ms/pair differ from incremental's.
* `--overwrite two_view_geometries` is the default: geometric verification
  does not rewrite raw matches, so both arms verify identical input without
  re-matching.
* Use a fresh `--run-name` when the configuration changes; a workspace caches
  the database it was built with, cameras included.
* colmap's output goes to
  `<run-dir>/eth3d/dslr/<scene>/{extraction,matching,reconstruction}.log`,
  which is also where each script's per-scene setup line is written.
