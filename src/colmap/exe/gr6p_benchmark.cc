// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Benchmark comparing optimized vs legacy GR6P estimator.

#include "colmap/estimators/generalized_relative_pose.h"
#include "colmap/estimators/generalized_relative_pose_legacy.h"
#include "colmap/estimators/poselib_utils.h"
#include "colmap/geometry/rigid3.h"
#include "colmap/math/random.h"
#include "colmap/optim/ransac.h"
#include "colmap/util/logging.h"
#include "colmap/util/timer.h"

#include <Eigen/Core>
#include <PoseLib/solvers/gen_relpose_6pt.h>

#include <iostream>
#include <vector>

using namespace colmap;

// Generate a synthetic relative pose problem.
void GenerateProblem(int num_points,
                     int num_cameras,
                     std::vector<GRNPObservation>* obs1,
                     std::vector<GRNPObservation>* obs2,
                     std::vector<GRNPObservationLegacy>* obs1_legacy,
                     std::vector<GRNPObservationLegacy>* obs2_legacy,
                     Rigid3d* rig2_from_rig1) {
  // Generate random rig poses
  const Rigid3d rig1_from_world(Eigen::Quaterniond::UnitRandom(),
                                Eigen::Vector3d::Random().normalized());
  const Rigid3d rig2_from_world(Eigen::Quaterniond::UnitRandom(),
                                Eigen::Vector3d::Random().normalized());
  *rig2_from_rig1 = rig2_from_world * Inverse(rig1_from_world);

  // Generate camera poses within each rig
  std::vector<Rigid3d> cams_from_rig1(num_cameras);
  std::vector<Rigid3d> cams_from_rig2(num_cameras);
  for (int i = 0; i < num_cameras; ++i) {
    cams_from_rig1[i] =
        Rigid3d(Eigen::Quaterniond::UnitRandom(),
                Eigen::Vector3d::Random().normalized());
    cams_from_rig2[i] =
        Rigid3d(Eigen::Quaterniond::UnitRandom(),
                Eigen::Vector3d::Random().normalized());
  }

  // Generate random 3D points and project them
  obs1->resize(num_points);
  obs2->resize(num_points);
  obs1_legacy->resize(num_points);
  obs2_legacy->resize(num_points);

  for (int i = 0; i < num_points; ++i) {
    const Eigen::Vector3d point3D = Eigen::Vector3d::Random() * 10;
    const int cam_idx1 = i % num_cameras;
    const int cam_idx2 = i % num_cameras;

    const Eigen::Vector3d point_in_cam1 =
        cams_from_rig1[cam_idx1] * (rig1_from_world * point3D);
    const Eigen::Vector3d point_in_cam2 =
        cams_from_rig2[cam_idx2] * (rig2_from_world * point3D);

    // Optimized observations
    (*obs1)[i].cam_from_rig = cams_from_rig1[cam_idx1].ToMatrix();
    (*obs1)[i].ray_in_cam = point_in_cam1.normalized();
    (*obs2)[i].cam_from_rig = cams_from_rig2[cam_idx2].ToMatrix();
    (*obs2)[i].ray_in_cam = point_in_cam2.normalized();

    // Legacy observations
    (*obs1_legacy)[i].cam_from_rig = cams_from_rig1[cam_idx1];
    (*obs1_legacy)[i].ray_in_cam = point_in_cam1.normalized();
    (*obs2_legacy)[i].cam_from_rig = cams_from_rig2[cam_idx2];
    (*obs2_legacy)[i].ray_in_cam = point_in_cam2.normalized();
  }
}

int main(int argc, char** argv) {
  InitializeGlog(argv);

  const int num_points = (argc >= 2) ? std::stoi(argv[1]) : 100;
  const int num_cameras = (argc >= 3) ? std::stoi(argv[2]) : 4;
  const int num_iterations = (argc >= 4) ? std::stoi(argv[3]) : 10000;

  std::cout << "Benchmark: GR6P Estimate Breakdown" << std::endl;
  std::cout << "  num_points: " << num_points << std::endl;
  std::cout << "  num_cameras: " << num_cameras << std::endl;
  std::cout << "  num_iterations: " << num_iterations << std::endl;

  SetPRNGSeed(42);

  // Generate a single problem
  std::vector<GRNPObservation> obs1, obs2;
  std::vector<GRNPObservationLegacy> obs1_legacy, obs2_legacy;
  Rigid3d rig2_from_rig1;

  GenerateProblem(num_points, num_cameras, &obs1, &obs2, &obs1_legacy,
                  &obs2_legacy, &rig2_from_rig1);

  // Convert obs1 to src format for optimized version
  std::vector<GRNPSrcObservation> src_obs1(num_points);
  for (int i = 0; i < num_points; ++i) {
    src_obs1[i] = obs1[i].ToSrc();
  }

  // Pick 6 random points for Estimate calls
  std::vector<GRNPSrcObservation> sample1(6);
  std::vector<GRNPObservation> sample2(6);
  std::vector<GRNPObservationLegacy> sample1_legacy(6);
  std::vector<GRNPObservationLegacy> sample2_legacy(6);
  for (int i = 0; i < 6; ++i) {
    sample1[i] = src_obs1[i];
    sample2[i] = obs2[i];
    sample1_legacy[i] = obs1_legacy[i];
    sample2_legacy[i] = obs2_legacy[i];
  }

  // Benchmark data prep only (optimized)
  double time_data_prep;
  {
    Timer timer;
    timer.Start();
    for (int t = 0; t < num_iterations; ++t) {
      std::vector<Eigen::Vector3d> origins_in_rig1(6);
      std::vector<Eigen::Vector3d> origins_in_rig2(6);
      std::vector<Eigen::Vector3d> rays_in_rig1(6);
      std::vector<Eigen::Vector3d> rays_in_rig2(6);
      for (int i = 0; i < 6; ++i) {
        const Eigen::Matrix3x4d& rfc1 = sample1[i].rig_from_cam;
        origins_in_rig1[i] = rfc1.col(3);
        rays_in_rig1[i] = rfc1.leftCols<3>() * sample1[i].ray_in_cam;

        const Eigen::Matrix3x4d& cfr2 = sample2[i].cam_from_rig;
        const Eigen::Matrix3d R2_inv = cfr2.leftCols<3>().transpose();
        origins_in_rig2[i] = -R2_inv * cfr2.col(3);
        rays_in_rig2[i] = R2_inv * sample2[i].ray_in_cam;
      }
    }
    timer.Pause();
    time_data_prep = timer.ElapsedSeconds();
  }

  // Benchmark poselib solver only
  double time_poselib;
  int total_poses = 0;
  {
    // Prepare data once
    std::vector<Eigen::Vector3d> origins_in_rig1(6);
    std::vector<Eigen::Vector3d> origins_in_rig2(6);
    std::vector<Eigen::Vector3d> rays_in_rig1(6);
    std::vector<Eigen::Vector3d> rays_in_rig2(6);
    for (int i = 0; i < 6; ++i) {
      const Eigen::Matrix3x4d& rfc1 = sample1[i].rig_from_cam;
      origins_in_rig1[i] = rfc1.col(3);
      rays_in_rig1[i] = rfc1.leftCols<3>() * sample1[i].ray_in_cam;

      const Eigen::Matrix3x4d& cfr2 = sample2[i].cam_from_rig;
      const Eigen::Matrix3d R2_inv = cfr2.leftCols<3>().transpose();
      origins_in_rig2[i] = -R2_inv * cfr2.col(3);
      rays_in_rig2[i] = R2_inv * sample2[i].ray_in_cam;
    }

    Timer timer;
    timer.Start();
    for (int t = 0; t < num_iterations; ++t) {
      std::vector<poselib::CameraPose> poses;
      poselib::gen_relpose_6pt(
          origins_in_rig1, rays_in_rig1, origins_in_rig2, rays_in_rig2, &poses);
      total_poses += poses.size();
    }
    timer.Pause();
    time_poselib = timer.ElapsedSeconds();
  }

  // Benchmark result conversion only
  double time_conversion;
  {
    // Get some poses first
    std::vector<Eigen::Vector3d> origins_in_rig1(6);
    std::vector<Eigen::Vector3d> origins_in_rig2(6);
    std::vector<Eigen::Vector3d> rays_in_rig1(6);
    std::vector<Eigen::Vector3d> rays_in_rig2(6);
    for (int i = 0; i < 6; ++i) {
      const Eigen::Matrix3x4d& rfc1 = sample1[i].rig_from_cam;
      origins_in_rig1[i] = rfc1.col(3);
      rays_in_rig1[i] = rfc1.leftCols<3>() * sample1[i].ray_in_cam;

      const Eigen::Matrix3x4d& cfr2 = sample2[i].cam_from_rig;
      const Eigen::Matrix3d R2_inv = cfr2.leftCols<3>().transpose();
      origins_in_rig2[i] = -R2_inv * cfr2.col(3);
      rays_in_rig2[i] = R2_inv * sample2[i].ray_in_cam;
    }
    std::vector<poselib::CameraPose> poses;
    poselib::gen_relpose_6pt(
        origins_in_rig1, rays_in_rig1, origins_in_rig2, rays_in_rig2, &poses);

    Timer timer;
    timer.Start();
    for (int t = 0; t < num_iterations; ++t) {
      std::vector<Rigid3d> results;
      results.reserve(poses.size());
      for (const poselib::CameraPose& pose : poses) {
        results.emplace_back(ConvertPoseLibPoseToRigid3d(pose));
      }
    }
    timer.Pause();
    time_conversion = timer.ElapsedSeconds();
  }

  // Benchmark full Estimate (optimized)
  double time_estimate_optimized;
  {
    Timer timer;
    timer.Start();
    for (int t = 0; t < num_iterations; ++t) {
      std::vector<Rigid3d> models;
      GR6PEstimator::Estimate(sample1, sample2, &models);
    }
    timer.Pause();
    time_estimate_optimized = timer.ElapsedSeconds();
  }

  // Benchmark full Estimate (legacy)
  double time_estimate_legacy;
  {
    Timer timer;
    timer.Start();
    for (int t = 0; t < num_iterations; ++t) {
      std::vector<Rigid3d> models;
      GR6PEstimatorLegacy::Estimate(sample1_legacy, sample2_legacy, &models);
    }
    timer.Pause();
    time_estimate_legacy = timer.ElapsedSeconds();
  }

  const double avg_poses = static_cast<double>(total_poses) / num_iterations;

  std::cout << "\nEstimate Breakdown (" << num_iterations << " iterations):"
            << std::endl;
  std::cout << "  Data prep:       " << (time_data_prep * 1000) << " ms ("
            << (time_data_prep / time_estimate_optimized * 100) << "%)"
            << std::endl;
  std::cout << "  Poselib solver:  " << (time_poselib * 1000) << " ms ("
            << (time_poselib / time_estimate_optimized * 100) << "%)"
            << std::endl;
  std::cout << "  Conversion:      " << (time_conversion * 1000) << " ms ("
            << (time_conversion / time_estimate_optimized * 100) << "%)"
            << std::endl;
  std::cout << "  Avg poses/call:  " << avg_poses << std::endl;

  std::cout << "\nFull Estimate:" << std::endl;
  std::cout << "  Optimized: " << (time_estimate_optimized * 1000) << " ms"
            << std::endl;
  std::cout << "  Legacy:    " << (time_estimate_legacy * 1000) << " ms"
            << std::endl;
  std::cout << "  Speedup:   " << (time_estimate_legacy / time_estimate_optimized)
            << "x" << std::endl;

  return EXIT_SUCCESS;
}
