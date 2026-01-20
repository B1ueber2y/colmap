// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Legacy GRNP estimator using Rigid3d for benchmarking comparison.

#pragma once

#include "colmap/geometry/rigid3.h"
#include "colmap/util/eigen_alignment.h"

#include <vector>

#include <Eigen/Core>

namespace colmap {

// Legacy observation struct using Rigid3d.
struct GRNPObservationLegacy {
  Rigid3d cam_from_rig;
  Eigen::Vector3d ray_in_cam;
};

// Legacy GR6P estimator using Rigid3d for benchmarking comparison.
class GR6PEstimatorLegacy {
 public:
  typedef GRNPObservationLegacy X_t;
  typedef GRNPObservationLegacy Y_t;
  typedef Rigid3d M_t;

  static const int kMinNumSamples = 6;

  static void Estimate(const std::vector<X_t>& points1,
                       const std::vector<Y_t>& points2,
                       std::vector<M_t>* rigs2_from_rigs1);

  static void Residuals(const std::vector<X_t>& points1,
                        const std::vector<Y_t>& points2,
                        const M_t& rig2_from_rig1,
                        std::vector<double>* residuals);
};

}  // namespace colmap
