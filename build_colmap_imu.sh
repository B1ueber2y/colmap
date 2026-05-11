#!/bin/bash
#SBATCH --job-name=build_colmap_imu
#SBATCH --partition=normal.24h
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --account=public
#SBATCH --time=4:00:00
#SBATCH --output=/cluster/scratch/ademirtas/spack_build_logs/%x_%j.out
#SBATCH --error=/cluster/scratch/ademirtas/spack_build_logs/%x_%j.err

set -e

module load eth_proxy
. /cluster/software/stacks/2025-06/setup-env.sh
source /cluster/scratch/ademirtas/code/colmap_imu/setup_env.sh

COLMAP_IMU_DIR=/cluster/scratch/ademirtas/code/colmap_imu
BUILD_DIR=$COLMAP_IMU_DIR/build
INSTALL_DIR=$COLMAP_IMU_DIR/install

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

CUDA_ROOT=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/cuda-12.6.2-6ovgc2ircbnyntakw37gaj2j7k2fuwp3
MESA_GLU=/cluster/software/stacks/2025-06/linux-ubuntu22.04-x86_64_v3/gcc-12.2.0/mesa-glu-9.0.2-64i4rzmaopawgynfbgh6xkgo4uctpo73

cmake .. -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCUDA_ENABLED=ON \
  -DCUDA_TOOLKIT_ROOT_DIR="$CUDA_ROOT" \
  -DCMAKE_CXX_FLAGS="-isystem $MESA_GLU/include" \
  -DCMAKE_CUDA_FLAGS="-isystem $MESA_GLU/include" \
  -DGUI_ENABLED=ON \
  -DOPENGL_ENABLED=ON \
  -DONNX_ENABLED=ON \
  -DCGAL_ENABLED=ON \
  -DTESTS_ENABLED=OFF \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR"

ninja -j $SLURM_CPUS_PER_TASK
ninja install

echo "Done: colmap_imu built at $BUILD_DIR/src/colmap/exe/colmap"
echo "      installed at $INSTALL_DIR/bin/colmap"
