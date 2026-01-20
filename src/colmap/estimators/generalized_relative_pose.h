// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "colmap/geometry/rigid3.h"
#include "colmap/util/eigen_alignment.h"

#include <vector>

#include <Eigen/Core>

namespace colmap {

// Forward declaration.
struct GRNPSrcObservation;

// General observation storing cam_from_rig as matrix.
struct GRNPObservation {
  Eigen::Matrix3x4d cam_from_rig;
  Eigen::Vector3d ray_in_cam;

  // Convert to source observation by computing the inverse.
  GRNPSrcObservation ToSrc() const;
};

// Source observation (points1) - stores rig_from_cam for fast Residuals.
struct GRNPSrcObservation {
  Eigen::Matrix3x4d rig_from_cam;
  Eigen::Vector3d ray_in_cam;
};

inline GRNPSrcObservation GRNPObservation::ToSrc() const {
  GRNPSrcObservation src;
  // Compute inverse: [R|t]^-1 = [R^T | -R^T*t]
  const Eigen::Matrix3d R_inv = cam_from_rig.leftCols<3>().transpose();
  src.rig_from_cam.leftCols<3>() = R_inv;
  src.rig_from_cam.col(3) = -R_inv * cam_from_rig.col(3);
  src.ray_in_cam = ray_in_cam;
  return src;
}

// Target observation (points2) - alias for GRNPObservation.
using GRNPTgtObservation = GRNPObservation;

// Minimal generalized relative pose estimator based on poselib.
class GR6PEstimator {
 public:
  typedef GRNPSrcObservation X_t;
  typedef GRNPTgtObservation Y_t;
  // The estimated rig2_from_rig1 relative pose between the generalized cameras.
  typedef Rigid3d M_t;

  // The minimum number of samples needed to estimate a model. Note that in
  // theory the minimum required number of samples is 6 but Laurent Kneip showed
  // in his paper that using 8 samples is more stable.
  static const int kMinNumSamples = 6;

  // Estimate the most probable solution of the GR6P problem from a set of
  // six 2D-2D point correspondences.
  static void Estimate(const std::vector<X_t>& points1,
                       const std::vector<Y_t>& points2,
                       std::vector<M_t>* rigs2_from_rigs1);

  // Calculate the squared Sampson error between corresponding points.
  static void Residuals(const std::vector<X_t>& points1,
                        const std::vector<Y_t>& points2,
                        const M_t& rig2_from_rig1,
                        std::vector<double>* residuals);
};

// Solver for the Generalized Relative Pose problem using a minimal of 8 2D-2D
// correspondences. This implementation is based on:
//
//    "Efficient Computation of Relative Pose for Multi-Camera Systems",
//    Kneip and Li. CVPR 2014.
//
// Note that the solution to this problem is degenerate in the case of pure
// translation and when all correspondences are observed from the same cameras.
//
// The implementation is a modified and improved version of Kneip's original
// implementation in OpenGV licensed under the BSD license.
class GR8PEstimator {
 public:
  typedef GRNPSrcObservation X_t;
  typedef GRNPTgtObservation Y_t;
  // The estimated rig2_from_rig1 relative pose between the generalized cameras.
  typedef Rigid3d M_t;

  // The minimum number of samples needed to estimate a model. Note that in
  // theory the minimum required number of samples is 6 but Laurent Kneip showed
  // in his paper that using 8 samples is more stable.
  static const int kMinNumSamples = 8;

  // Estimate the most probable solution of the GR6P problem from a set of
  // six 2D-2D point correspondences.
  static void Estimate(const std::vector<X_t>& points1,
                       const std::vector<Y_t>& points2,
                       std::vector<M_t>* rigs2_from_rigs1);

  // Calculate the squared Sampson error between corresponding points.
  static void Residuals(const std::vector<X_t>& points1,
                        const std::vector<Y_t>& points2,
                        const M_t& rig2_from_rig1,
                        std::vector<double>* residuals);
};

}  // namespace colmap
