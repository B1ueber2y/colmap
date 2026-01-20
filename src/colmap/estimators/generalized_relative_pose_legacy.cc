// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Legacy GRNP estimator using Rigid3d for benchmarking comparison.

#include "colmap/estimators/generalized_relative_pose_legacy.h"

#include "colmap/estimators/poselib_utils.h"
#include "colmap/geometry/essential_matrix.h"
#include "colmap/util/logging.h"

#include <PoseLib/solvers/gen_relpose_6pt.h>

namespace colmap {

void GR6PEstimatorLegacy::Estimate(const std::vector<X_t>& points1,
                                   const std::vector<Y_t>& points2,
                                   std::vector<M_t>* rigs2_from_rigs1) {
  THROW_CHECK_EQ(points1.size(), 6);
  THROW_CHECK_EQ(points2.size(), 6);
  THROW_CHECK_NOTNULL(rigs2_from_rigs1);

  rigs2_from_rigs1->clear();

  std::vector<Eigen::Vector3d> origins_in_rig1(6);
  std::vector<Eigen::Vector3d> origins_in_rig2(6);
  std::vector<Eigen::Vector3d> rays_in_rig1(6);
  std::vector<Eigen::Vector3d> rays_in_rig2(6);
  for (int i = 0; i < 6; ++i) {
    origins_in_rig1[i] = points1[i].cam_from_rig.TgtOriginInSrc();
    origins_in_rig2[i] = points2[i].cam_from_rig.TgtOriginInSrc();
    rays_in_rig1[i] =
        points1[i].cam_from_rig.rotation.inverse() * points1[i].ray_in_cam;
    rays_in_rig2[i] =
        points2[i].cam_from_rig.rotation.inverse() * points2[i].ray_in_cam;
  }

  std::vector<poselib::CameraPose> poses;
  poselib::gen_relpose_6pt(
      origins_in_rig1, rays_in_rig1, origins_in_rig2, rays_in_rig2, &poses);

  rigs2_from_rigs1->reserve(poses.size());
  for (const poselib::CameraPose& pose : poses) {
    rigs2_from_rigs1->emplace_back(ConvertPoseLibPoseToRigid3d(pose));
  }
}

void GR6PEstimatorLegacy::Residuals(const std::vector<X_t>& points1,
                                    const std::vector<Y_t>& points2,
                                    const M_t& rig2_from_rig1,
                                    std::vector<double>* residuals) {
  THROW_CHECK_EQ(points1.size(), points2.size());
  residuals->resize(points1.size());
  for (size_t i = 0; i < points1.size(); ++i) {
    const Rigid3d cam2_from_cam1 = points2[i].cam_from_rig * rig2_from_rig1 *
                                   Inverse(points1[i].cam_from_rig);
    const Eigen::Matrix3d E = EssentialMatrixFromPose(cam2_from_cam1);
    const Eigen::Vector3d epipolar_line1 = E * points1[i].ray_in_cam;
    const double num = points2[i].ray_in_cam.dot(epipolar_line1);
    const Eigen::Vector4d denom(points2[i].ray_in_cam.dot(E.col(0)),
                                points2[i].ray_in_cam.dot(E.col(1)),
                                epipolar_line1.x(),
                                epipolar_line1.y());
    (*residuals)[i] = num * num / denom.squaredNorm();
  }
}

}  // namespace colmap
