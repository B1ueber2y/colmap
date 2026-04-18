# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

COLMAP is a Structure-from-Motion (SfM) and Multi-View Stereo (MVS) pipeline for 3D reconstruction from image collections. This repository is on the `features/imu` branch, which adds IMU (Inertial Measurement Unit) integration for Visual-Inertial Odometry (VIO).

## Build Commands

```bash
# Standard build
mkdir build && cd build
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install
ninja

# Minimal build (no GUI, no CUDA — faster for development)
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release -DGUI_ENABLED=OFF -DCUDA_ENABLED=OFF

# Build with tests
cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release -DTESTS_ENABLED=ON

# Python bindings (after C++ install)
colmap_DIR=./install ./python/incremental_build.sh   # incremental
colmap_DIR=./install ./python/build.sh               # clean
```

## Running Tests

```bash
cd build
ctest --output-on-failure             # all tests
ctest -R "sensor/imu_test"            # specific test
ctest -R "estimators/imu"             # pattern match
```

Python tests: `pytest` from repo root.

Test files are named `*_test.cc` and registered via `COLMAP_ADD_TEST()`. CTest names follow `module/test_name` (e.g., `estimators/alignment_test`).

## Formatting & Linting

```bash
scripts/format/c++.sh     # clang-format (Google style) + clang-tidy
scripts/format/python.sh  # ruff (line length 80)
```

CI enforces both. Fix format issues before pushing — CI runs clang-format and ruff checks as separate jobs.

## Architecture

### Module Dependency Layers (bottom → top)

| Layer | Modules |
|-------|---------|
| Primitives | `util/`, `math/` |
| Sensor | `sensor/` (camera models, IMU calibration), `geometry/` |
| Data | `scene/` (Camera, Image, Frame, Point3D, Reconstruction, ImuState) |
| Features | `feature/` (SIFT, ALIKED, LightGlue) |
| Optimization | `optim/` (RANSAC), `estimators/` (BA, pose, IMU preintegration) |
| Pipelines | `sfm/`, `mvs/`, `controllers/` |
| Interface | `exe/` (CLI), `ui/` (Qt GUI), `pycolmap/` (Python bindings) |

### Key Data Structures

- **`Reconstruction`** — top-level container holding cameras, images, frames, 3D points, tracks
- **`Camera`** — intrinsics + distortion model (SimplePinhole, Radial, OpenCV, Fisheye, etc.)
- **`Image`** — 2D observations + pose, linked to a `Frame`
- **`Frame`** — posed rig instantiation; multiple images can share a frame
- **`Database`** — SQLite store for features and matches
- **`ImuState`** — 9-DOF state: velocity (3) + gyro bias (3) + accel bias (3)

### IMU Integration (features/imu)

The IMU additions span four modules:

**`src/colmap/sensor/imu.h/.cc`** — `ImuCalibration` (noise/bias params), `ImuMeasurement`, `ImuMeasurements` (sorted, thread-safe). Based on ADIS16448 defaults.

**`src/colmap/scene/imu.h`** — `Imu` struct (calibration + `imu_from_cam` transform), `ImuState` (contiguous param block with Eigen::Map accessors).

**`src/colmap/estimators/imu_preintegration.h/.cc`** — `ImuPreintegrator` supporting two methods:
- `MIDPOINT`: trapezoidal, fast, for high-rate IMUs
- `RK4` (default): analytical rotation integrals + RK4 covariance propagation

Produces `PreintegratedImuData` with delta R/p/v, bias Jacobians, and a 15×15 covariance (condition-number clamped `sqrt_information`). `ImuReintegrationCallback` integrates with Ceres to reintegrate when bias drift exceeds threshold.

**`src/colmap/estimators/cost_functions/imu_preintegration.h`** — Two Ceres cost functions with 15D residuals (rotation error, position error, velocity error, bias random walk):
- `ImuPreintegrationCostFunctor`: AutoDiff
- `AnalyticalImuPreintegrationCostFunction`: hand-derived Jacobians (Forster et al. TRO 2016)

Parameter blocks: `body_from_world_i` (7), `imu_state_i` (9), `body_from_world_j` (7), `imu_state_j` (9).

**`python/examples/vi_optimization.py`** — Full VIO optimization example showing calibration setup, measurement loading, preintegration, and Ceres BA with IMU factors.

### Naming Conventions

- Classes/methods: `PascalCase`
- Member variables: `snake_case_` (trailing underscore)
- Transforms: `target_from_source` (e.g., `imu_from_cam`, `body_from_world_i`)
- Coordinates: `point_in_frame` (e.g., `point3D_in_world`)
- Special ID types: `camera_t`, `image_t`, `frame_t`, `point3D_t`, `sensor_t` (from `util/types.h`)

### Adding New Modules

Use the `COLMAP_ADD_LIBRARY()` and `COLMAP_ADD_TEST()` CMake macros defined in `cmake/CMakeHelper.cmake`. Tests go in `*_test.cc` files alongside source files.
