// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.

#include "colmap/estimators/generalized_pose.h"
#include "colmap/geometry/rigid3.h"
#include "colmap/optim/ransac.h"
#include "colmap/scene/camera.h"
#include "colmap/sensor/models.h"
#include "colmap/util/logging.h"
#include "colmap/util/timer.h"

#include <Eigen/Core>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace colmap;

int main(int argc, char** argv) {
  InitializeGlog(argv);

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <rig_test_data.txt> [max_trials]"
              << std::endl;
    return EXIT_FAILURE;
  }

  const std::string input_path = argv[1];
  const int max_trials = (argc >= 3) ? std::stoi(argv[2]) : 1000;

  std::ifstream file(input_path);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << input_path << std::endl;
    return EXIT_FAILURE;
  }

  std::string line;
  auto read_line = [&]() {
    while (std::getline(file, line)) {
      if (!line.empty() && line[0] != '#') {
        return true;
      }
    }
    return false;
  };

  // Read threshold
  read_line();
  const double thresh = std::stod(line);
  std::cout << "Threshold: " << thresh << std::endl;

  // Read number of cameras
  read_line();
  const size_t num_cameras = std::stoul(line);
  std::cout << "Number of cameras: " << num_cameras << std::endl;

  // Read cameras
  std::vector<Camera> cameras(num_cameras);
  for (size_t i = 0; i < num_cameras; ++i) {
    read_line();
    std::istringstream iss(line);
    std::string model_name;
    size_t width, height;
    iss >> model_name >> width >> height;

    cameras[i].model_id = CameraModelNameToId(model_name);
    cameras[i].width = width;
    cameras[i].height = height;
    cameras[i].params.resize(CameraModelNumParams(cameras[i].model_id));
    for (size_t j = 0; j < cameras[i].params.size(); ++j) {
      iss >> cameras[i].params[j];
    }
    std::cout << "Camera " << i << ": " << model_name << " " << width << "x"
              << height << std::endl;
  }

  // Read cams_from_rig transforms
  std::vector<Rigid3d> cams_from_rig(num_cameras);
  for (size_t i = 0; i < num_cameras; ++i) {
    read_line();
    std::istringstream iss(line);
    double qw, qx, qy, qz, tx, ty, tz;
    iss >> qw >> qx >> qy >> qz >> tx >> ty >> tz;
    // Eigen quaternion uses xyzw internally
    cams_from_rig[i].rotation = Eigen::Quaterniond(qw, qx, qy, qz);
    cams_from_rig[i].translation = Eigen::Vector3d(tx, ty, tz);
  }

  // Read number of correspondences
  read_line();
  const size_t num_correspondences = std::stoul(line);
  std::cout << "Number of correspondences: " << num_correspondences << std::endl;

  // Read correspondences
  std::vector<Eigen::Vector2d> points2D(num_correspondences);
  std::vector<Eigen::Vector3d> points3D(num_correspondences);
  std::vector<size_t> camera_idxs(num_correspondences);

  for (size_t i = 0; i < num_correspondences; ++i) {
    read_line();
    std::istringstream iss(line);
    double x2d, y2d, x3d, y3d, z3d;
    iss >> camera_idxs[i] >> x2d >> y2d >> x3d >> y3d >> z3d;
    points2D[i] = Eigen::Vector2d(x2d, y2d);
    points3D[i] = Eigen::Vector3d(x3d, y3d, z3d);
  }

  file.close();

  // Setup RANSAC options
  RANSACOptions ransac_options;
  ransac_options.max_error = thresh;
  ransac_options.max_num_trials = max_trials;
  ransac_options.min_num_trials = max_trials;
  ransac_options.confidence = 0.9999;

  std::cout << "\nRunning EstimateGeneralizedAbsolutePose with max_trials="
            << max_trials << std::endl;

  // Run estimation
  Rigid3d rig_from_world;
  size_t num_inliers;
  std::vector<char> inlier_mask;

  Timer timer;
  timer.Start();
  const bool success = EstimateGeneralizedAbsolutePose(ransac_options,
                                                       points2D,
                                                       points3D,
                                                       camera_idxs,
                                                       cams_from_rig,
                                                       cameras,
                                                       &rig_from_world,
                                                       &num_inliers,
                                                       &inlier_mask);
  timer.Pause();

  std::cout << "\nResults:" << std::endl;
  std::cout << "  Success: " << (success ? "true" : "false") << std::endl;
  std::cout << "  Num inliers: " << num_inliers << std::endl;
  std::cout << "  Total time: " << timer.ElapsedSeconds() << " s" << std::endl;
  std::cout << "  Time per trial: " << (timer.ElapsedSeconds() * 1000 / max_trials)
            << " ms" << std::endl;

  if (success) {
    std::cout << "  Rotation (qw qx qy qz): " << rig_from_world.rotation.w()
              << " " << rig_from_world.rotation.x() << " "
              << rig_from_world.rotation.y() << " "
              << rig_from_world.rotation.z() << std::endl;
    std::cout << "  Translation: " << rig_from_world.translation.transpose()
              << std::endl;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
