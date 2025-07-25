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

#include "colmap/util/eigen_alignment.h"
#include "colmap/util/types.h"

#include <vector>

#include <Eigen/Core>
#include <ceres/ceres.h>

namespace colmap {

// Essential matrix estimator from corresponding normalized camera ray pairs.
//
// This algorithm solves the 5-Point problem based on the following paper:
//
//    D. Nister, An efficient solution to the five-point relative pose problem,
//    IEEE-T-PAMI, 26(6), 2004.
//    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.86.8769
class EssentialMatrixFivePointEstimator {
 public:
  typedef Eigen::Vector3d X_t;
  typedef Eigen::Vector3d Y_t;
  typedef Eigen::Matrix3d M_t;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 5;

  // Estimate up to 10 possible essential matrix solutions from a set of
  // corresponding camera rays.
  //
  //  The number of corresponding rays must be at least 5.
  //
  // @param cam_rays1  First set of corresponding rays.
  // @param cam_rays2  Second set of corresponding rays.
  //
  // @return           Up to 10 solutions as a vector of 3x3 essential matrices.
  static void Estimate(const std::vector<X_t>& cam_rays1,
                       const std::vector<Y_t>& cam_rays2,
                       std::vector<M_t>* models);

  // Calculate the residuals of a set of corresponding rays and a given
  // essential matrix.
  //
  // Residuals are defined as the squared Sampson error.
  //
  // @param cam_rays1  First set of corresponding rays.
  // @param cam_rays2  Second set of corresponding rays.
  // @param E          3x3 essential matrix.
  // @param residuals  Output vector of residuals.
  static void Residuals(const std::vector<X_t>& cam_rays1,
                        const std::vector<Y_t>& cam_rays2,
                        const M_t& E,
                        std::vector<double>* residuals);
};

// Essential matrix estimator from corresponding normalized camera ray pairs.
//
// This algorithm solves the 8-Point problem based on the following paper:
//
//    Hartley and Zisserman, Multiple View Geometry, algorithm 11.1, page 282.
class EssentialMatrixEightPointEstimator {
 public:
  typedef Eigen::Vector3d X_t;
  typedef Eigen::Vector3d Y_t;
  typedef Eigen::Matrix3d M_t;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 8;

  // Estimate essential matrix solutions from set of corresponding camera rays.

  //
  // The number of corresponding rays must be at least 8.
  //
  // @param cam_rays1  First set of corresponding rays.
  // @param cam_rays2  Second set of corresponding rays.
  static void Estimate(const std::vector<X_t>& cam_rays1,
                       const std::vector<Y_t>& cam_rays2,

                       std::vector<M_t>* models);

  // Calculate the residuals of a set of corresponding rays and a given
  // essential matrix.
  //
  // Residuals are defined as the squared Sampson error.
  //
  // @param cam_rays1  First set of corresponding rays.
  // @param cam_rays2  Second set of corresponding rays.
  // @param E          3x3 essential matrix.
  // @param residuals  Output vector of residuals.
  static void Residuals(const std::vector<X_t>& cam_rays1,
                        const std::vector<Y_t>& cam_rays2,
                        const M_t& E,
                        std::vector<double>* residuals);
};

}  // namespace colmap
