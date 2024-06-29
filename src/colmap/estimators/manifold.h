// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
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

#include <ceres/ceres.h>
#include <ceres/rotation.h>

namespace colmap {

#include <cmath>

#include "ceres/ceres.h"

inline void SetQuaternionManifold(ceres::Problem* problem, double* quat_xyzw) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
  problem->SetManifold(quat_xyzw, new ceres::EigenQuaternionManifold);
#else
  problem->SetParameterization(quat_xyzw,
                               new ceres::EigenQuaternionParameterization);
#endif
}

inline void SetSubsetManifold(int size,
                              const std::vector<int>& constant_params,
                              ceres::Problem* problem,
                              double* params) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
  problem->SetManifold(params,
                       new ceres::SubsetManifold(size, constant_params));
#else
  problem->SetParameterization(
      params, new ceres::SubsetParameterization(size, constant_params));
#endif
}

template <int size>
inline void SetSphereManifold(ceres::Problem* problem, double* params) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
  problem->SetManifold(params, new ceres::SphereManifold<size>);
#else
  problem->SetParameterization(
      params, new ceres::HomogeneousVectorParameterization(size));
#endif
}

// Use an exponential function to ensure the variable to be strictly positive
// Generally for the scale parameters
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
template <int size = 1>
class PositiveExponentialManifold : public ceres::Manifold {
 public:
  PositiveExponentialManifold() {}
  ~PositiveExponentialManifold() {}

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override {
    for (size_t i = 0; i < size; ++i) {
      x_plus_delta[i] = std::exp(std::log(x[i]) + delta[i]);
    }
    return true;
  }

  bool PlusJacobian(const double* x, double* jacobian) const override {
    for (size_t i = 0; i < size; ++i) {
      jacobian[size * i + i] = x[i];
    }
    return true;
  }

  int AmbientSize() const override { return size; }
  int TangentSize() const override { return size; }
};
#else
template <int size = 1>
class PositiveExponentialParameterization
    : public ceres::LocalParameterization {
 public:
  PositiveExponentialParameterization() {}
  ~PositiveExponentialParameterization() {}

  bool Plus(const double* x,
            const double* delta,
            double* x_plus_delta) const override {
    for (size_t i = 0; i < size; ++i) {
      x_plus_delta[i] = std::exp(std::log(x[i]) + delta[i]);
    }
    return true;
  }

  bool ComputeJacobian(const double* x, double* jacobian) const override {
    for (size_t i = 0; i < size; ++i) {
      jacobian[size * i + i] = x[i];
    }
    return true;
  }

  int GlobalSize() const override { return size; }
  int LocalSize() const override { return size; }
};
#endif

template <int size = 1>
inline void SetPositiveExponentialManifold(ceres::Problem* problem,
                                           double* params) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
  problem->SetManifold(params, new PositiveExponentialManifold<size>);
#else
  problem->SetParameterization(quat_xyzw,
                               new PositiveExponentialParameterization<size>);
#endif
}

inline int ParameterBlockTangentSize(ceres::Problem* problem,
                                     const double* param) {
#if CERES_VERSION_MAJOR >= 3 || \
    (CERES_VERSION_MAJOR == 2 && CERES_VERSION_MINOR >= 1)
  return problem->ParameterBlockTangentSize(param);
#else
  return problem->ParameterBlockLocalSize(param);
#endif
}

}  // namespace colmap
