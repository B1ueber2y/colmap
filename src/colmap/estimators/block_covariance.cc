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

#include "colmap/estimators/block_covariance.h"

#include "colmap/estimators/covariance.h"

namespace colmap {

ProblemPartitioner::ProblemPartitioner(ceres::Problem* problem,
                                       Reconstruction* reconstruction) {
  problem_ = problem;
  for (const auto& image : reconstruction->Images()) {
    const double* qvec = image.second.CamFromWorld().rotation.coeffs().data();
    if (problem_->HasParameterBlock(qvec))
      pose_blocks_.insert(const_cast<double*>(qvec));
    const double* tvec = image.second.CamFromWorld().translation.data();
    if (problem_->HasParameterBlock(tvec))
      pose_blocks_.insert(const_cast<double*>(tvec));
  }
  BuildBipartiteGraph();
}

ProblemPartitioner::ProblemPartitioner(
    ceres::Problem* problem, const std::vector<const double*>& pose_blocks) {
  problem_ = problem;
  for (auto it = pose_blocks.begin(); it != pose_blocks.end(); ++it) {
    if (problem->HasParameterBlock(*it))
      pose_blocks_.insert(const_cast<double*>(*it));
  }
  BuildBipartiteGraph();
}

void ProblemPartitioner::BuildBipartiteGraph() {
  std::vector<ceres::ResidualBlockId> residual_block_ids;
  problem_->GetResidualBlocks(&residual_block_ids);
  for (const auto& residual_block_id : residual_block_ids) {
    std::vector<double*> param_blocks;
    problem_->GetParameterBlocksForResidualBlock(residual_block_id,
                                                 &param_blocks);

    for (auto& param_block : param_blocks) {
      graph_.AddEdge(param_block, residual_block_id);
    }
  }
}

void ProblemPartitioner::GetSubproblem(
    ceres::Problem* subproblem,
    const std::vector<double*>& subproblem_pose_blocks) const {
  THROW_CHECK_EQ(subproblem->NumResiduals(), 0);
  // construct irrelevant pose blocks
  std::unordered_set<double*> irrelevant_pose_blocks;
  for (auto it = pose_blocks_.begin(); it != pose_blocks_.end(); ++it) {
    if (std::find(subproblem_pose_blocks.begin(),
                  subproblem_pose_blocks.end(),
                  *it) != subproblem_pose_blocks.end())
      continue;
    irrelevant_pose_blocks.insert(*it);
  }

  // customizd bfs to collect relevant residuals
  std::unordered_set<const double*> param_visited;
  std::unordered_set<ceres::ResidualBlockId> residual_block_ids;
  for (const auto& param : subproblem_pose_blocks) {
    THROW_CHECK(problem_->HasParameterBlock(param));
    if (param_visited.find(param) != param_visited.end()) continue;
    std::queue<double*> bfs_queue;
    bfs_queue.push(param);
    param_visited.insert(param);
    while (!bfs_queue.empty()) {
      double* current_param = bfs_queue.front();
      bfs_queue.pop();
      for (auto& residual_block_id : graph_.GetResidualBlocks(current_param)) {
        // check if the residual block exists
        if (residual_block_ids.find(residual_block_id) !=
            residual_block_ids.end())
          continue;
        // skip if the residual block links to irrelevant poses
        bool flag = true;
        for (auto& param_block : graph_.GetParameterBlocks(residual_block_id)) {
          if (irrelevant_pose_blocks.find(param_block) !=
              irrelevant_pose_blocks.end()) {
            flag = false;
            break;
          }
        }
        if (!flag) continue;
        // insert residual blocks
        residual_block_ids.insert(residual_block_id);
        for (auto& param_block : graph_.GetParameterBlocks(residual_block_id)) {
          if (param_visited.find(param_block) != param_visited.end()) continue;
          bfs_queue.push(param_block);
          param_visited.insert(param_block);
        }
      }
    }
  }

  // build subproblem
  for (auto& residual_block_id : residual_block_ids) {
    subproblem->AddResidualBlock(
        problem_->residual_block(residual_block_id)->release());
  }
}

LocationPartitioner::LocationPartitioner(
    const std::map<image_t, Eigen::Vector2d>& locations)
    : locations_(locations) {
  std::vector<BGPointWithIndex> points;
  for (auto it = locations.begin(); it != locations.end(); ++it) {
    points.emplace_back(BGPoint(it->second(0), it->second(1)), it->first);
  }
  rtree_.insert(points.begin(), points.end());
}

LocationPartitioner::LocationPartitioner(Reconstruction* reconstruction) {
  std::vector<BGPointWithIndex> points;
  for (const auto& image : reconstruction->Images()) {
    Eigen::Vector3d center = image.second.ProjectionCenter();
    locations_.emplace(image.first, Eigen::Vector2d(center(0), center(1)));
    points.emplace_back(BGPoint(center(0), center(1)), image.first);
  }
  rtree_.insert(points.begin(), points.end());
}

std::vector<image_t> LocationPartitioner::GetImageIdsInsideRange(
    const Eigen::Vector2d& top_left,
    const Eigen::Vector2d& bottom_right) const {
  THROW_CHECK_LT(top_left(0), bottom_right(0));
  THROW_CHECK_LT(top_left(1), bottom_right(1));

  BGPoint bgpoint_top_left(top_left(0), top_left(1));
  BGPoint bgpoint_bottom_right(bottom_right(0), bottom_right(1));
  boost::geometry::model::box<BGPoint> query_box(bgpoint_top_left,
                                                 bgpoint_bottom_right);
  std::vector<BGPointWithIndex> result;
  rtree_.query(boost::geometry::index::intersects(query_box),
               std::back_inserter(result));

  std::vector<image_t> output;
  for (auto it = result.begin(); it != result.end(); ++it) {
    output.push_back(it->second);
  }
  return output;
}

std::vector<double> LocationPartitioner::GetCoords(
    const std::vector<double>& inputs,
    double block_size,
    double robust_percentile,
    double kStretchRatio) const {
  // sort list
  std::vector<double> xs;
  xs.insert(xs.end(), inputs.begin(), inputs.end());
  std::sort(xs.begin(), xs.end());
  size_t n_images = xs.size();

  // compute the robust ranges (x_start and x_end)
  double robust_x_min = xs[int(robust_percentile * n_images)];
  double robust_x_max = xs[int((1 - robust_percentile) * n_images)];
  double center = 0.5 * (robust_x_min + robust_x_max);
  int n_coords = int(std::ceil((robust_x_max - robust_x_min) *
                               (1. + kStretchRatio) / block_size));
  double x_start = center - n_coords * block_size / 2.0;
  double x_end = center + n_coords * block_size / 2.0;

  // construct coords
  std::vector<double> coords;
  if (xs[0] < x_start) coords.push_back(xs[0]);
  for (int i = 0; i < n_coords; ++i) {
    coords.push_back(x_start + block_size * i);
  }
  if (xs[n_images - 1] > x_end) coords.push_back(xs[n_images - 1]);
  return coords;
}

std::vector<LocationPartitioner::Box2d> LocationPartitioner::GetPartitions(
    double block_size, double robust_percentile, double kStretchRatio) const {
  std::vector<double> xs, ys;
  for (auto it = locations_.begin(); it != locations_.end(); ++it) {
    xs.push_back(it->second(0));
    ys.push_back(it->second(1));
  }
  std::vector<double> coords_x =
      GetCoords(xs, block_size, robust_percentile, kStretchRatio);
  std::vector<double> coords_y =
      GetCoords(ys, block_size, robust_percentile, kStretchRatio);
  std::vector<Box2d> boxes;
  for (size_t i = 0; i < coords_x.size() - 1; ++i) {
    for (size_t j = 0; j < coords_y.size() - 1; ++j) {
      Eigen::Vector2d top_left(coords_x[i], coords_y[j]);
      Eigen::Vector2d bottom_right(coords_x[i + 1], coords_y[j + 1]);
      boxes.push_back(std::make_pair(top_left, bottom_right));
    }
  }
  return boxes;
}

bool BlockwiseCovarianceEstimator::Estimate(
    std::map<image_t, Eigen::MatrixXd>& image_id_to_covar) {
  ProblemPartitioner problem_partitioner(problem_, reconstruction_);
  LocationPartitioner location_partitioner(reconstruction_);

  std::vector<LocationPartitioner::Box2d> boxes =
      location_partitioner.GetPartitions(options_.block_size,
                                         options_.robust_percentile,
                                         options_.kStretchRatio);

  double margin = (options_.window_size - options_.block_size) / 2.0;
  for (const auto& box2d : boxes) {
    Eigen::Vector2d top_left = box2d.first;
    Eigen::Vector2d bottom_right = box2d.second;
    std::vector<image_t> image_ids =
        location_partitioner.GetImageIdsInsideRange(top_left, bottom_right);
    // stretch with margin
    Eigen::Vector2d window_top_left =
        top_left - Eigen::Vector2d(margin, margin);
    Eigen::Vector2d window_bottom_right =
        bottom_right + Eigen::Vector2d(margin, margin);
    std::vector<image_t> window_image_ids =
        location_partitioner.GetImageIdsInsideRange(top_left, bottom_right);
    // get subproblem
    std::vector<double*> subproblem_pose_blocks;
    for (const auto& image : reconstruction_->Images()) {
      const double* qvec = image.second.CamFromWorld().rotation.coeffs().data();
      if (problem_->HasParameterBlock(qvec))
        subproblem_pose_blocks.push_back(const_cast<double*>(qvec));
      const double* tvec = image.second.CamFromWorld().translation.data();
      if (problem_->HasParameterBlock(tvec))
        subproblem_pose_blocks.push_back(const_cast<double*>(tvec));
    }
    ceres::Problem* subproblem = new ceres::Problem();
    problem_partitioner.GetSubproblem(subproblem, subproblem_pose_blocks);
    std::map<image_t, Eigen::MatrixXd> tmp_image_id_to_covar;
    if (!EstimatePoseCovariance(
            subproblem, reconstruction_, tmp_image_id_to_covar))
      return false;
    for (const image_t& image_id : image_ids) {
      if (image_id_to_covar.find(image_id) != image_id_to_covar.end()) continue;
      image_id_to_covar.emplace(image_id, tmp_image_id_to_covar.at(image_id));
    }
  }
  return true;
}

}  // namespace colmap
