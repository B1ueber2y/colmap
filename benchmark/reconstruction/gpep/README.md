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

All numbers below are **`--mapper incremental`**. The global pipeline is not
measured here; see the note on it under Notes.

Measured on one RTX 3080 at `--num-parallel-scenes 1`, five seeds, both arms
from one worktree at `af20f0fa` differing only in `kEnableFocalSolvers`.
`nofocal` is the without-arm, `focal` is stock main, and `A - B` is
`focal - nofocal`, so positive favours the solver.

### Runtime

| incremental mapper | verification, without | verification, with | change |
| --- | --- | --- | --- |
| shared focal | 9.9 ms/pair | 95.7 ms/pair | +862% |
| one-sided focal | 18.7 ms/pair | 46.8 ms/pair | +150% |

One seed, reusing the accuracy run's workspaces. The denominator is the pairs
verification ran on (those holding raw matches): 10662 and 15422. The figure
is the matching phase total, which at `--quality high` also includes guided
matching.

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
* `--mapper global` works but measures a different thing. The incremental
  mapper consumes the focal these solvers estimate
  (`incremental_mapper_impl.cc`, `info.camera1 = two_view_geometry.camera1`),
  whereas view-graph calibration recomputes focals from the view graph and
  then clears `tvg.camera1/camera2` so consumers use its own K. Under the
  global mapper the solvers therefore only affect the pipeline indirectly,
  through the F they produce and the inliers they keep. The global path also
  disables guided matching and uses different two-view thresholds
  (`max_error` 1.0, `min_num_inliers` 30, `min_inlier_ratio` 0.25), so its
  runtime figure isolates the solver far better than the incremental one.
* `--overwrite two_view_geometries` is the default: geometric verification
  does not rewrite raw matches, so both arms verify identical input without
  re-matching.
* Use a fresh `--run-name` when the configuration changes; a workspace caches
  the database it was built with, cameras included.
* colmap's output goes to
  `<run-dir>/eth3d/dslr/<scene>/{extraction,matching,reconstruction}.log`,
  which is also where each script's per-scene setup line is written.
