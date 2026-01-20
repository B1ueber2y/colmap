// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Legacy GP3P estimator using Rigid3d for benchmarking comparison.

#include "colmap/estimators/generalized_absolute_pose_legacy.h"

#include "colmap/estimators/poselib_utils.h"
#include "colmap/util/logging.h"

#include <PoseLib/solvers/gp3p.h>
#include <PoseLib/solvers/p3p.h>

namespace colmap {

GP3PEstimatorLegacy::GP3PEstimatorLegacy(ResidualType residual_type)
    : residual_type_(residual_type) {}

void GP3PEstimatorLegacy::Estimate(const std::vector<X_t>& points2D,
                                   const std::vector<Y_t>& points3D,
                                   std::vector<M_t>* rigs_from_world) const {
  THROW_CHECK_EQ(points2D.size(), 3);
  THROW_CHECK_EQ(points3D.size(), 3);
  THROW_CHECK_NOTNULL(rigs_from_world);

  rigs_from_world->clear();

  std::vector<Eigen::Vector3d> rays_in_rig(3);
  std::vector<Eigen::Vector3d> origins_in_rig(3);
  for (int i = 0; i < 3; ++i) {
    const Rigid3d rig_from_cam = Inverse(points2D[i].cam_from_rig);
    rays_in_rig[i] =
        (rig_from_cam.rotation * points2D[i].ray_in_cam).normalized();
    origins_in_rig[i] = rig_from_cam.translation;
  }

  std::vector<poselib::CameraPose> poses;
  if (origins_in_rig[0].isApprox(origins_in_rig[1], 1e-6) &&
      origins_in_rig[0].isApprox(origins_in_rig[2], 1e-6)) {
    // In case of a panoramic camera/rig, fall back to P3P.
    poselib::p3p(rays_in_rig, points3D, &poses);
    for (poselib::CameraPose& pose : poses) {
      pose.t += origins_in_rig[0];
    }
  } else {
    poselib::gp3p(origins_in_rig, rays_in_rig, points3D, &poses);
  }

  rigs_from_world->reserve(poses.size());
  for (const poselib::CameraPose& pose : poses) {
    rigs_from_world->emplace_back(ConvertPoseLibPoseToRigid3d(pose));
  }
}

void GP3PEstimatorLegacy::Residuals(const std::vector<X_t>& points2D,
                                    const std::vector<Y_t>& points3D,
                                    const M_t& rig_from_world,
                                    std::vector<double>* residuals) const {
  THROW_CHECK_EQ(points2D.size(), points3D.size());
  residuals->resize(points2D.size(), 0);
  for (size_t i = 0; i < points2D.size(); ++i) {
    const Eigen::Vector3d point3D_in_cam =
        points2D[i].cam_from_rig * (rig_from_world * points3D[i]);
    // Check if 3D point is in front of camera.
    if (point3D_in_cam.z() > std::numeric_limits<double>::epsilon()) {
      if (residual_type_ == ResidualType::CosineDistance) {
        const double cosine_dist =
            1 - point3D_in_cam.normalized().dot(points2D[i].ray_in_cam);
        (*residuals)[i] = cosine_dist * cosine_dist;
      } else if (residual_type_ == ResidualType::ReprojectionError) {
        (*residuals)[i] = (point3D_in_cam.hnormalized() -
                           points2D[i].ray_in_cam.hnormalized())
                              .squaredNorm();
      } else {
        LOG(FATAL_THROW) << "Invalid residual type";
      }
    } else {
      (*residuals)[i] = std::numeric_limits<double>::max();
    }
  }
}

}  // namespace colmap
