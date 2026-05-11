# `colmap_imu` Fork — Capabilities and Implementation Guide

**Date**: 2026-04-22
**Repo**: `/cluster/scratch/ademirtas/code/colmap_imu` (branch: `features/imu`)
**Context**: Assessment of which improvements from `COLMAP_PIPELINE_PROPOSALS.md` are directly
implementable using this fork, and how to implement them.

---

## Summary

The fork adds two independent capability sets on top of standard COLMAP:

1. **Post-hoc Visual-Inertial Refinement** — fully implemented, ready to use. Adds metric
   scale, gravity alignment, velocity, and bias estimation to an existing SfM reconstruction
   using IMU preintegration factors. This is the highest-value addition.

2. **ALIKED + LightGlue feature extraction/matching** — fully implemented, gated on the
   `COLMAP_ONNX_ENABLED` build flag and ONNX model files. Replaces SIFT for better performance
   on low-texture, repetitive-structure environments like construction sites.

IMU-guided preprocessing (frame selection, initial pair seeding, spatial matching pairs) is
**not** implemented in the fork but can be built as Python scripts using the exposed
`ImuPreintegrator` bindings.

---

## Key Files

| File | Purpose |
|---|---|
| `src/colmap/sensor/imu.h/.cc` | `ImuCalibration`, `ImuMeasurement`, `ImuMeasurements` |
| `src/colmap/scene/imu.h` | `Imu` struct, `ImuState` (velocity + biases) |
| `src/colmap/estimators/imu_preintegration.h/.cc` | `ImuPreintegrator` (MIDPOINT / RK4) |
| `src/colmap/estimators/cost_functions/imu_preintegration.h` | 4 Ceres cost functions |
| `src/colmap/feature/aliked.h` | ALIKED extractor + ALIKED+LightGlue matcher |
| `src/colmap/feature/onnx_matchers.h` | `LightGlueONNXFeatureMatcher`, `BruteForceONNXFeatureMatcher` |
| `src/pycolmap/estimators/imu_preintegration.cc` | Python bindings for all IMU classes |
| `src/pycolmap/feature/matching.cc` | Python bindings: `FeatureMatcherType` enum |
| `python/examples/vi_optimization.py` | Complete end-to-end VI refinement example |

---

## Capability 1: Post-Hoc Visual-Inertial Refinement

### What it does

Takes a finished COLMAP reconstruction (the `sparse/0_localized_refined/` output) and
refines it by adding IMU preintegration factors between consecutive registered frames.
Jointly optimizes:

| Variable | Dimensionality | Effect |
|---|---|---|
| `log_scale` | 1 | Recovers metric scale from IMU accelerometer |
| `gravity` | 3 (unit vector) | Aligns reconstruction to physical gravity direction |
| `imu_from_cam` | 7 (rigid transform) | Refines IMU-camera extrinsic calibration |
| `imu_state[i].velocity` | 3 per frame | Per-frame velocity, enforces physical continuity |
| `imu_state[i].bias_gyro` | 3 per frame | Gyroscope bias |
| `imu_state[i].bias_accel` | 3 per frame | Accelerometer bias |

The cost function used is `AnalyticalVisualCentricImuPreintegrationCostFunction` (hand-derived
Jacobians, Forster et al. TRO 2016), which is designed for starting from a visual-only SfM in
an arbitrary frame — exactly the situation after running the current COLMAP pipeline.

### Why this matters for the Hilti challenge

The current pipeline already has metric scale from the stereo rig baseline, but:
- Gravity is not explicitly constrained — the reconstruction sits in an arbitrary orientation
- Velocity consistency across frames is not enforced — pure SfM has no dynamics model
- The Hilti GT is gravity-aligned; any residual tilt in the COLMAP frame contributes to ATE

VI refinement forces the trajectory to be physically consistent: an agent with these poses
must have plausible accelerations and angular rates given the IMU readings. This is an
orthogonal constraint to reprojection that can correct slow drift not visible in BA cost.

### The rig blocker and workaround

**Blocker**: The VI cost functions assert `len(image.frame.rig.non_ref_sensors) == 0`. The
current pipeline uses a stereo rig (cam0 + cam1 as a rigid pair), which fails this assert.
The IMU is not yet integrated into the `Rig` abstraction (explicit TODO in `scene/imu.h:47`).

**Workaround (no C++ changes required)**:

Extract cam0-only poses from the reconstruction, run VI refinement on those, then recompose
the stereo poses using the fixed rig transform:

```python
import pycolmap
import numpy as np

# 1. Load the full stereo reconstruction
recon = pycolmap.Reconstruction("sparse_txt/0_localized_refined/")

# 2. Build a cam0-only sub-reconstruction (trivial rig)
recon_cam0 = pycolmap.Reconstruction()
for image_id, image in recon.images.items():
    if "cam0" in image.name:
        recon_cam0.add_image(image)
        # ... (copy cameras, frames, points as needed)

# 3. Run VI refinement on recon_cam0
# (see vi_optimization.py for the full setup)
run_vi_optimization(recon_cam0, database_path, output_path, integrators, preintegrated, variables)

# 4. Recompose cam1 poses using the known fixed rig transform
# cam1_from_world = cam1_from_cam0 * cam0_from_world
rig_transform = load_rig_config("rig_config.json")  # known from factory calibration
for image in refined_recon.images.values():
    if "cam0" in image.name:
        cam1_image = find_corresponding_cam1(image)
        cam1_image.cam_from_world = rig_transform * image.cam_from_world
```

This loses nothing: cam1 poses are fully determined by cam0 + the fixed rig transform, which
is already how the current pipeline treats them (rig extrinsics frozen).

### Hilti IMU calibration

The Hilti sensor uses a **Bosch BMI085** IMU. Calibration values to use (replace the Project
Aria defaults in `vi_optimization.py`):

```python
imu_calib = pycolmap.ImuCalibration()
imu_calib.imu_rate = 200.0                        # Hz
imu_calib.gravity_magnitude = 9.81
imu_calib.gyro_noise_density = 1.7e-4             # rad/s/sqrt(Hz) — BMI085 spec
imu_calib.accel_noise_density = 2.0e-3            # m/s^2/sqrt(Hz) — BMI085 spec
imu_calib.bias_gyro_random_walk_sigma = 1.9e-5    # rad/s^2/sqrt(Hz)
imu_calib.bias_accel_random_walk_sigma = 3.0e-4   # m/s^3/sqrt(Hz)
```

Verify against the Hilti dataset's `imu_params.yaml` if available; these are nominal datasheet
values.

### Timestamp extraction from Insta360 filenames

Filenames follow `XXXXXX_NNNNNNNNNNNNNNN.png` where `N` is nanoseconds. Extract as:

```python
def filename_to_ns(name: str) -> int:
    return int(Path(name).stem.split("_")[1])
```

### Adding VI refinement as Stage 9 in the pipeline

```bash
# After stage 8 (model converter), add:
python run_vi_refinement.py \
    --reconstruction  cfg_10hz_jointBA_robust/sparse/0_localized_refined/ \
    --database        cfg_10hz_jointBA_robust/database.db \
    --imu_data        /path/to/hilti/imu_measurements.npy \
    --image_timestamps /path/to/image_timestamps.npy \
    --output          cfg_10hz_vi_refined/sparse/
```

---

## Capability 2: ALIKED + LightGlue Features

### Build requirement

ALIKED and LightGlue are compiled only when `COLMAP_ONNX_ENABLED=ON`:

```bash
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release \
    -DONNX_ENABLED=ON \
    -DONNXRUNTIME_INCLUDE_DIR=/path/to/onnxruntime/include \
    -DONNXRUNTIME_LIBRARY=/path/to/onnxruntime/lib/libonnxruntime.so
```

ONNX model files are downloaded from default URIs at runtime
(`kDefaultAlikedN16RotFeatureExtractorUri`, `kDefaultAlikedLightGlueFeatureMatcherUri`).
Check `src/colmap/feature/resources.h` for the exact URLs.

### Matcher types available (from Python bindings)

```python
import pycolmap

# Four matcher types exposed:
pycolmap.FeatureMatcherType.SIFT_BRUTEFORCE     # standard COLMAP
pycolmap.FeatureMatcherType.SIFT_LIGHTGLUE      # ← best quick win: keeps SIFT features
pycolmap.FeatureMatcherType.ALIKED_BRUTEFORCE   # ALIKED + brute-force cosine similarity
pycolmap.FeatureMatcherType.ALIKED_LIGHTGLUE    # ← best quality: full learned pipeline
```

### Recommended option for the Hilti pipeline

**`SIFT_LIGHTGLUE` is the quickest win**: it replaces the mutual-NN SIFT matching with
LightGlue while keeping existing SIFT feature extractions. This means:
- No need to re-extract features for all 290K images (expensive)
- Just re-run matching stage with LightGlue
- LightGlue's cross-attention reduces outlier ratio from ~90% (SIFT NN on construction)
  to ~20–30%, which means `abs_pose_min_inlier_ratio=0.10` can likely be raised back
  toward 0.25 and `abs_pose_min_num_inliers=15` can be raised toward 20+

**`ALIKED_LIGHTGLUE` is the ceiling**: ALIKED is more robust than SIFT on textureless
surfaces (concrete, plywood shuttering, bare metal). Requires full re-extraction — plan
~2× the current feature extraction time per run.

### Using LightGlue in the sequential matching stage

The existing SLURM scripts call `colmap sequential_matcher`. To swap in LightGlue, the
`FeatureMatcherType` needs to be set. Check whether the CLI exposes this via
`--SequentialMatching.matcher_type` or whether it must be done through pycolmap's
sequential matching controller.

```python
# Via pycolmap (if CLI flag not available):
matching_options = pycolmap.SequentialMatchingOptions()
matching_options.matcher_type = pycolmap.FeatureMatcherType.SIFT_LIGHTGLUE
pycolmap.match_sequential(database_path, matching_options)
```

### LightGlue options

```python
# From src/colmap/feature/onnx_matchers.h
lightglue_opts = pycolmap.LightGlueONNXMatchingOptions()
lightglue_opts.min_score = 0.1   # post-model confidence threshold
# lower = more matches (but noisier), higher = fewer but cleaner

# For construction site (high outlier ratio with SIFT keypoints):
lightglue_opts.min_score = 0.15  # slightly stricter than default
```

---

## Capability 3: IMU Preprocessing Scripts (Not in Fork — Build on Top)

These three scripts do not exist yet but can be written in ~1–2 days each using the
`pycolmap.ImuPreintegrator` bindings that are already there.

### Script A: Motion-adaptive frame selection

Replaces the fixed `10hz_frames.txt` (every 3rd frame) with a velocity-gated list.

**Input**: Hilti IMU data (200 Hz), image timestamps  
**Output**: `adaptive_frames.txt` (same format as `10hz_frames.txt`)

```python
import pycolmap
import numpy as np
from pathlib import Path
from scipy.spatial import KDTree

def build_adaptive_frame_list(
    imu_measurements: pycolmap.ImuMeasurements,
    image_timestamps: dict[int, int],  # frame_idx → nanoseconds
    imu_calib: pycolmap.ImuCalibration,
    min_translation_m: float = 0.05,
    max_rotation_deg: float = 20.0,
    max_angular_vel_rad_s: float = 1.5,
    output_path: Path = Path("adaptive_frames.txt"),
) -> list[int]:
    opts = pycolmap.ImuPreintegrationOptions()
    selected = []
    last_selected_ts = None

    sorted_frames = sorted(image_timestamps.items())
    for i, (frame_idx, ts) in enumerate(sorted_frames):
        if last_selected_ts is None:
            selected.append(frame_idx)
            last_selected_ts = ts
            continue

        ms = imu_measurements.extract_measurements_contain_edge(last_selected_ts, ts)
        if len(ms) == 0:
            continue

        integrator = pycolmap.ImuPreintegrator(opts, imu_calib, last_selected_ts, ts)
        integrator.feed_imu(ms)
        data = integrator.extract()

        translation_norm = np.linalg.norm(data.delta_p)
        rotation_deg = np.degrees(2 * np.arccos(
            min(1.0, abs(data.delta_R.w))
        ))

        # Estimate angular velocity from last measurement in window
        last_gyro = ms[-1].gyro if hasattr(ms[-1], 'gyro') else np.zeros(3)
        angular_vel = np.linalg.norm(last_gyro)

        if (angular_vel < max_angular_vel_rad_s
                and translation_norm > min_translation_m
                and rotation_deg < max_rotation_deg):
            selected.append(frame_idx)
            last_selected_ts = ts

    output_path.write_text("\n".join(str(f) for f in selected))
    return selected
```

### Script B: IMU-seeded initial pair selection

Finds frame pairs with good triangulation geometry, passes them to the COLMAP mapper CLI.

**Input**: Preintegrated IMU poses  
**Output**: ranked list of `(image_id1, image_id2)` candidates

```python
def find_initial_pair_candidates(
    imu_measurements: pycolmap.ImuMeasurements,
    image_timestamps: dict[int, int],
    imu_calib: pycolmap.ImuCalibration,
    min_baseline_m: float = 0.3,
    max_baseline_m: float = 2.0,
    max_rotation_deg: float = 30.0,
    top_k: int = 5,
) -> list[tuple[int, int]]:
    opts = pycolmap.ImuPreintegrationOptions()
    frame_list = sorted(image_timestamps.items())  # [(frame_idx, ts), ...]
    candidates = []

    for i in range(len(frame_list)):
        for j in range(i + 1, min(i + 50, len(frame_list))):
            t1, t2 = frame_list[i][1], frame_list[j][1]
            ms = imu_measurements.extract_measurements_contain_edge(t1, t2)
            if len(ms) == 0:
                continue
            integrator = pycolmap.ImuPreintegrator(opts, imu_calib, t1, t2)
            integrator.feed_imu(ms)
            data = integrator.extract()

            baseline = np.linalg.norm(data.delta_p)
            rot_deg = np.degrees(2 * np.arccos(min(1.0, abs(data.delta_R.w))))

            if min_baseline_m <= baseline <= max_baseline_m and rot_deg < max_rotation_deg:
                score = baseline / (rot_deg + 1e-3)  # prefer translation over rotation
                candidates.append((score, frame_list[i][0], frame_list[j][0]))

    candidates.sort(reverse=True)
    return [(a, b) for _, a, b in candidates[:top_k]]

# Then pass to COLMAP:
# best_pair = candidates[0]
# colmap mapper ... --Mapper.init_image_id1 {best_pair[0]} --Mapper.init_image_id2 {best_pair[1]}
```

### Script C: Spatial matching pairs from IMU poses

Accumulates full-trajectory positions via preintegration, builds a k-d tree, finds spatially
close but temporally distant pairs for `colmap custom_matcher`.

**Input**: Full-trajectory preintegrated positions  
**Output**: `spatial_pairs.txt` (image_name_i image_name_j per line)

```python
def build_spatial_pairs(
    imu_measurements: pycolmap.ImuMeasurements,
    image_timestamps: dict[int, int],
    image_names: dict[int, str],
    imu_calib: pycolmap.ImuCalibration,
    max_spatial_dist_m: float = 2.0,
    min_temporal_gap_frames: int = 100,
    output_path: Path = Path("spatial_pairs.txt"),
) -> None:
    opts = pycolmap.ImuPreintegrationOptions()
    frame_list = sorted(image_timestamps.items())

    # Accumulate positions via sequential preintegration
    positions = np.zeros((len(frame_list), 3))
    for i in range(1, len(frame_list)):
        t1, t2 = frame_list[i-1][1], frame_list[i][1]
        ms = imu_measurements.extract_measurements_contain_edge(t1, t2)
        if len(ms) == 0:
            positions[i] = positions[i-1]
            continue
        integrator = pycolmap.ImuPreintegrator(opts, imu_calib, t1, t2)
        integrator.feed_imu(ms)
        data = integrator.extract()
        positions[i] = positions[i-1] + data.delta_p

    # Spatial pairs via k-d tree
    tree = KDTree(positions)
    pairs = set()
    for i, (frame_idx, _) in enumerate(frame_list):
        nearby = tree.query_ball_point(positions[i], r=max_spatial_dist_m)
        for j in nearby:
            if abs(i - j) >= min_temporal_gap_frames:
                key = (min(i, j), max(i, j))
                pairs.add(key)

    with open(output_path, "w") as f:
        for i, j in sorted(pairs):
            f.write(f"{image_names[frame_list[i][0]]} {image_names[frame_list[j][0]]}\n")

# Then run:
# colmap custom_matcher --database_path database.db \
#     --CustomMatching.match_list_path spatial_pairs.txt \
#     --CustomMatching.match_type pairs
```

---

## Integration into the Existing Pipeline

### New Stage 0.5: IMU Preprocessing (before feature extraction)

```bash
# Runs once per sequence, produces 3 files:
python imu_preprocessing.py \
    --imu_data      /path/to/hilti_imu.npy \
    --timestamps    /path/to/image_timestamps.npy \
    --image_dir     data_30Hz/floor_X/run_Y/images/rig1/cam0/ \
    --output_dir    ./ \
    # Writes: adaptive_frames.txt, init_pair_candidates.txt, spatial_pairs.txt
```

### Modified Stage 2: Matching (spatial pairs added)

```bash
# After sequential matcher:
colmap custom_matcher \
    --database_path database.db \
    --CustomMatching.match_list_path spatial_pairs.txt \
    --CustomMatching.match_type pairs
```

### Modified Stage 3: Mapper (seeded initialization)

```bash
# Read top candidate from init_pair_candidates.txt and pass:
colmap mapper \
    ... \
    --Mapper.init_image_id1 NNNN \
    --Mapper.init_image_id2 MMMM \
    # (other existing flags unchanged)
```

### New Stage 9: Visual-Inertial Refinement (after existing Stage 8)

```bash
python stage9_vi_refinement.py \
    --reconstruction cfg_10hz_jointBA_robust/sparse/0_localized_refined/ \
    --database       cfg_10hz_jointBA_robust/database.db \
    --imu_data       /path/to/hilti_imu.npy \
    --timestamps     /path/to/image_timestamps.npy \
    --output         cfg_10hz_vi_refined/sparse/
# Then run model converter on cfg_10hz_vi_refined/sparse/
```

---

## What Requires C++ Changes (Deferred)

| Feature | Blocker | Effort |
|---|---|---|
| VI refinement on stereo rig directly | `scene/imu.h:47` TODO — IMU not in Rig abstraction | ~3 days C++ |
| IMU factors during incremental SfM | Not wired into `IncrementalMapper` | ~1 week C++ |
| IMU-guided sequential matcher | Matcher has no pose prior input | ~2 days C++ |

The cam0-only VI refinement workaround unblocks the first item without C++ changes.

---

## Recommended Implementation Order

1. **VI refinement (Stage 9)**: highest impact, Python-only, reference implementation already
   exists in `vi_optimization.py`. Adapt for Hilti IMU format and cam0-only rig workaround.
   Test on one completed GT run to measure scale/gravity correction magnitude.

2. **SIFT + LightGlue matching**: if `COLMAP_ONNX_ENABLED` is on in the build, swap
   `FeatureMatcherType` to `SIFT_LIGHTGLUE` in the sequential matching stage. No re-extraction.
   Run on one sequence, compare coverage and BA cost.

3. **Spatial matching pairs (Script C)**: ~1 day. Directly addresses the loop-closure gap for
   construction revisits. Low risk — additional pairs can only help, never hurt coverage.

4. **Adaptive frame selection (Script A)**: ~1 day. Replaces `10hz_frames.txt` generation.
   Compare coverage and BA cost; fall back to time-based if quality drops.

5. **IMU initial pair seeding (Script B)**: ~half day. Eliminates the init-pair lottery on
   hard sequences. Low risk — pass candidates via CLI flag, mapper still verifies geometry.
