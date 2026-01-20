// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Legacy GP3P estimator using Rigid3d for benchmarking comparison.

#pragma once

#include "colmap/geometry/rigid3.h"
#include "colmap/util/eigen_alignment.h"

#include <vector>

#include <Eigen/Core>

namespace colmap {

// Legacy solver for the Generalized P3P problem using Rigid3d.
// This is used for benchmarking against the optimized Matrix3x4d version.
class GP3PEstimatorLegacy {
 public:
  // The generalized image observations, which is composed of the relative pose
  // of a camera in the generalized camera and a ray in the camera frame.
  struct X_t {
    Rigid3d cam_from_rig;
    Eigen::Vector3d ray_in_cam;
  };

  // The observed 3D feature points in the world frame.
  typedef Eigen::Vector3d Y_t;
  // The estimated rig_from_world pose of the generalized camera.
  typedef Rigid3d M_t;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 3;

  // Whether to compute the cosine similarity or the reprojection error.
  enum class ResidualType {
    CosineDistance,
    ReprojectionError,
  };

  explicit GP3PEstimatorLegacy(
      ResidualType residual_type = ResidualType::CosineDistance);

  // Estimate the most probable solution of the GP3P problem from a set of
  // three 2D-3D point correspondences.
  void Estimate(const std::vector<X_t>& points2D,
                const std::vector<Y_t>& points3D,
                std::vector<M_t>* models) const;

  // Calculate the squared cosine distance error between the rays given a set of
  // 2D-3D point correspondences and the rig pose of the generalized camera.
  void Residuals(const std::vector<X_t>& points2D,
                 const std::vector<Y_t>& points3D,
                 const M_t& rig_from_world,
                 std::vector<double>* residuals) const;

 private:
  const ResidualType residual_type_;
};

}  // namespace colmap
