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

#include "colmap/scene/reconstruction.h"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <ceres/ceres.h>

namespace colmap {

class ProblemPartitioner {
 public:
  ProblemPartitioner(ceres::Problem* problem, Reconstruction* reconstruction);
  ProblemPartitioner(ceres::Problem* problem,
                     const std::vector<const double*>& pose_blocks);

  void GetSubproblem(ceres::Problem* subproblem,
                     const std::vector<double*>& subproblem_pose_blocks) const;

 private:
  struct BipartiteGraph {
    BipartiteGraph() {}

    void AddEdge(double* param_block,
                 ceres::ResidualBlockId residual_block_id) {
      param_to_residual[param_block].push_back(residual_block_id);
      residual_to_param[residual_block_id].push_back(param_block);
    }
    std::vector<ceres::ResidualBlockId> GetResidualBlocks(
        double* param_block) const {
      return param_to_residual.at(param_block);
    }
    std::vector<double*> GetParameterBlocks(
        ceres::ResidualBlockId residual_block_id) const {
      return residual_to_param.at(residual_block_id);
    }
    std::unordered_map<double*, std::vector<ceres::ResidualBlockId>>
        param_to_residual;
    std::unordered_map<ceres::ResidualBlockId, std::vector<double*>>
        residual_to_param;
  };

  ceres::Problem* problem_;
  std::set<double*> pose_blocks_;

  // for ceres problem
  void BuildBipartiteGraph();
  BipartiteGraph graph_;
};

class LocationPartitioner {
 public:
  LocationPartitioner(const std::map<image_t, Eigen::Vector2d>& locations);
  LocationPartitioner(Reconstruction* reconstruction);

  std::vector<image_t> GetImageIdsInsideRange(
      const Eigen::Vector2d& top_left,
      const Eigen::Vector2d& bottom_right) const;

  using Box2d =
      std::pair<Eigen::Vector2d, Eigen::Vector2d>;  // (top_left, bottom_right)
  std::vector<Box2d> GetPartitions(double block_size,
                                   double robust_percentile = 0.01,
                                   double kStretchRatio = 0.1) const;

 private:
  // locations
  std::map<image_t, Eigen::Vector2d> locations_;

  // RTree data structure
  typedef boost::geometry::model::
      point<double, 2, boost::geometry::cs::cartesian>
          BGPoint;
  typedef std::pair<BGPoint, image_t> BGPointWithIndex;
  boost::geometry::index::rtree<BGPointWithIndex,
                                boost::geometry::index::quadratic<16>>
      rtree_;

  // get coords on 1D
  std::vector<double> GetCoords(const std::vector<double>& inputs,
                                double block_size,
                                double robust_percentile = 0.01,
                                double kStretchRatio = 0.1) const;
};

struct BlockwiseCovarianceEstimatorOptions {
  // The 2D window size for each subproblem
  double window_size = 300.0;
  // The 2D block size (inside the window) for all interested parameters in each
  // subproblem
  double block_size = 100.0;

  // Robust ranges and stretching
  double robust_percentile = 0.01;
  double kStretchRatio = 0.1;
};

class BlockwiseCovarianceEstimator {
 public:
  // Construct with a COLMAP reconstruction
  BlockwiseCovarianceEstimator(BlockwiseCovarianceEstimatorOptions& options,
                               ceres::Problem* problem,
                               Reconstruction* reconstruction,
                               double lambda = 1e-8)
      : options_(options),
        problem_(problem),
        reconstruction_(reconstruction),
        lambda_(lambda) {}

  // Run blockwise estimation
  bool Estimate(std::map<image_t, Eigen::MatrixXd>& image_id_to_covar);

 private:
  // Options
  const BlockwiseCovarianceEstimatorOptions options_;

  // ceres problem
  ceres::Problem* problem_;

  // reconstruction
  Reconstruction* reconstruction_;

  // The damping factor to avoid rank deficiency
  const double lambda_ = 1e-8;
};

}  // namespace colmap
