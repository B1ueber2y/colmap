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

#include "colmap/math/math.h"
#include "colmap/util/eigen_alignment.h"
#include "colmap/util/enum_utils.h"
#include "colmap/util/types.h"

#include <array>
#include <cfloat>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/LU>
#include <ceres/jet.h>

namespace colmap {

// This file defines several different camera models and arbitrary new camera
// models can be added by the following steps:
//
//  1. Add a new struct in this file which implements all the necessary methods.
//  2. Define an unique model_name and model_id for the camera model.
//  3. Add camera model to `CAMERA_MODEL_CASES` macro in this file.
//  4. Add new template specialization of test case for camera model to
//     `camera_models_test.cc`.
//
// A camera model can have three different types of camera parameters: focal
// length, principal point, extra parameters (distortion parameters). The
// parameter array is split into different groups, so that we can enable or
// disable the refinement of the individual groups during bundle adjustment. It
// is up to the camera model to access the parameters correctly (it is free to
// do so in an arbitrary manner) - the parameters are not accessed from outside.
//
// A camera model must have the following methods:
//
//  - `ImgFromCam`: Projects points in camera frame to pixel coordinates in
//    image plane (the inverse of `CamFromImg`). Assumes that the camera
//    coordinates are given as (u, v, w). Returns false, if projection failed.
//  - `CamFromImg`: lift pixel coordinates in image plane to normalized camera
//    coordinates (the inverse of `ImgFromCam`). Produces camera coordinates
//    as (u, v, 1). Returns false, if lifting failed.
//  - `CamFromImgThreshold`: transform a threshold given in pixels to
//    normalized units (e.g. useful for reprojection error thresholds).
//
// Whenever you specify the camera parameters in a list, they must appear
// exactly in the order as they are accessed in the defined model struct.
//
// The camera models follow the convention that the upper left image corner has
// the coordinate (0, 0), the lower right corner (width, height), i.e. that
// the upper left pixel center has coordinate (0.5, 0.5) and the lower right
// pixel center has the coordinate (width - 0.5, height - 0.5).

MAKE_ENUM_CLASS_OVERLOAD_STREAM(CameraModelId,
                                -1,
                                kInvalid,                // = -1
                                kSimplePinhole,          // = 0
                                kPinhole,                // = 1
                                kSimpleRadial,           // = 2
                                kRadial,                 // = 3
                                kOpenCV,                 // = 4
                                kOpenCVFisheye,          // = 5
                                kFullOpenCV,             // = 6
                                kFOV,                    // = 7
                                kSimpleRadialFisheye,    // = 8
                                kRadialFisheye,          // = 9
                                kThinPrismFisheye,       // = 10
                                kRadTanThinPrismFisheye  // = 11
);

#ifndef CAMERA_MODEL_DEFINITIONS
#define CAMERA_MODEL_DEFINITIONS(model_id_val,                                \
                                 model_name_val,                              \
                                 num_focal_params_val,                        \
                                 num_pp_params_val,                           \
                                 num_extra_params_val)                        \
  static constexpr size_t num_params =                                        \
      (num_focal_params_val) + (num_pp_params_val) + (num_extra_params_val);  \
  static constexpr size_t num_focal_params = num_focal_params_val;            \
  static constexpr size_t num_pp_params = num_pp_params_val;                  \
  static constexpr size_t num_extra_params = num_extra_params_val;            \
  static constexpr CameraModelId model_id = model_id_val;                     \
  static const std::string model_name;                                        \
  static const std::string params_info;                                       \
  static const std::array<size_t, (num_focal_params_val)> focal_length_idxs;  \
  static const std::array<size_t, (num_pp_params_val)> principal_point_idxs;  \
  static const std::array<size_t, (num_extra_params_val)> extra_params_idxs;  \
                                                                              \
  static inline CameraModelId InitializeModelId() { return model_id_val; };   \
  static inline std::string InitializeModelName() { return model_name_val; }; \
  static inline std::string InitializeParamsInfo();                           \
  static inline std::array<size_t, (num_focal_params_val)>                    \
  InitializeFocalLengthIdxs();                                                \
  static inline std::array<size_t, (num_pp_params_val)>                       \
  InitializePrincipalPointIdxs();                                             \
  static inline std::array<size_t, (num_extra_params_val)>                    \
  InitializeExtraParamsIdxs();                                                \
                                                                              \
  static inline std::vector<double> InitializeParams(                         \
      double focal_length, size_t width, size_t height);                      \
  template <typename T>                                                       \
  static bool ImgFromCam(const T* params, T u, T v, T w, T* x, T* y);         \
  static inline bool CamFromImg(                                              \
      const double* params, double x, double y, double* u, double* v);        \
  template <typename T>                                                       \
  static void Distortion(const T* extra_params, T u, T v, T* du, T* dv);
#endif

#ifndef CAMERA_MODEL_CASES
#define CAMERA_MODEL_CASES                          \
  CAMERA_MODEL_CASE(SimplePinholeCameraModel)       \
  CAMERA_MODEL_CASE(PinholeCameraModel)             \
  CAMERA_MODEL_CASE(SimpleRadialCameraModel)        \
  CAMERA_MODEL_CASE(SimpleRadialFisheyeCameraModel) \
  CAMERA_MODEL_CASE(RadialCameraModel)              \
  CAMERA_MODEL_CASE(RadialFisheyeCameraModel)       \
  CAMERA_MODEL_CASE(OpenCVCameraModel)              \
  CAMERA_MODEL_CASE(OpenCVFisheyeCameraModel)       \
  CAMERA_MODEL_CASE(FullOpenCVCameraModel)          \
  CAMERA_MODEL_CASE(FOVCameraModel)                 \
  CAMERA_MODEL_CASE(ThinPrismFisheyeCameraModel)    \
  CAMERA_MODEL_CASE(RadTanThinPrismFisheyeModel)
#endif

#ifndef CAMERA_MODEL_SWITCH_CASES
#define CAMERA_MODEL_SWITCH_CASES         \
  CAMERA_MODEL_CASES                      \
  default:                                \
    CAMERA_MODEL_DOES_NOT_EXIST_EXCEPTION \
    break;
#endif

#define CAMERA_MODEL_DOES_NOT_EXIST_EXCEPTION \
  throw std::domain_error("Camera model does not exist");

// Fisheye camera model macros
#ifndef FISHEYE_CAMERA_MODEL_CASES
#define FISHEYE_CAMERA_MODEL_CASES                  \
  CAMERA_MODEL_CASE(SimpleRadialFisheyeCameraModel) \
  CAMERA_MODEL_CASE(RadialFisheyeCameraModel)       \
  CAMERA_MODEL_CASE(OpenCVFisheyeCameraModel)       \
  CAMERA_MODEL_CASE(ThinPrismFisheyeCameraModel)    \
  CAMERA_MODEL_CASE(RadTanThinPrismFisheyeModel)
#endif

#ifndef FISHEYE_CAMERA_MODEL_DEFINITIONS
#define FISHEYE_CAMERA_MODEL_DEFINITIONS                               \
  template <typename T>                                                \
  static void ImgFromFisheye(const T* params, T uu, T vv, T* x, T* y); \
  template <typename T>                                                \
  static void FisheyeFromImg(const T* params, T x, T y, T* uu, T* vv);
#endif

// The "Curiously Recurring Template Pattern" (CRTP) is used here, so that we
// can reuse some shared functionality between all camera models -
// defined in the BaseCameraModel.
template <typename CameraModel>
struct BaseCameraModel {
  template <typename T>
  static inline bool HasBogusParams(const std::vector<T>& params,
                                    size_t width,
                                    size_t height,
                                    T min_focal_length_ratio,
                                    T max_focal_length_ratio,
                                    T max_extra_param);

  template <typename T>
  static inline bool HasBogusFocalLength(const std::vector<T>& params,
                                         size_t width,
                                         size_t height,
                                         T min_focal_length_ratio,
                                         T max_focal_length_ratio);

  template <typename T>
  static inline bool HasBogusPrincipalPoint(const std::vector<T>& params,
                                            size_t width,
                                            size_t height);

  template <typename T>
  static inline bool HasBogusExtraParams(const std::vector<T>& params,
                                         T max_extra_param);

  template <typename T>
  static inline T CamFromImgThreshold(const T* params, T threshold);

  static inline bool IterativeUndistortion(const double* params,
                                           double* u,
                                           double* v);
};

// Base model for Fisheye camera models
template <typename CameraModel>
struct BaseFisheyeCameraModel : public BaseCameraModel<CameraModel> {
  template <typename T>
  static inline void FisheyeFromNormal(const T u, const T v, T* uu, T* vv) {
    *uu = u;
    *vv = v;
    const T r = ceres::sqrt(u * u + v * v);
    if (r > T(std::numeric_limits<double>::epsilon())) {
      const T theta = ceres::atan(r);
      *uu *= theta / r;
      *vv *= theta / r;
    }
  }

  template <typename T>
  static inline void NormalFromFisheye(const T uu, const T vv, T* u, T* v) {
    *u = uu;
    *v = vv;
    const T theta = ceres::sqrt(uu * uu + vv * vv);
    const T theta_cos_theta = theta * ceres::cos(theta);
    if (theta_cos_theta > T(std::numeric_limits<double>::epsilon())) {
      const T scale = ceres::sin(theta) / theta_cos_theta;
      *u *= scale;
      *v *= scale;
    }
  }
};

// Simple Pinhole camera model.
//
// No Distortion is assumed. Only focal length and principal point is modeled.
//
// Parameter list is expected in the following order:
//
//   f, cx, cy
//
// See https://en.wikipedia.org/wiki/Pinhole_camera_model
struct SimplePinholeCameraModel
    : public BaseCameraModel<SimplePinholeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kSimplePinhole, "SIMPLE_PINHOLE", 1, 2, 0)
};

// Pinhole camera model.
//
// No Distortion is assumed. Only focal length and principal point is modeled.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy
//
// See https://en.wikipedia.org/wiki/Pinhole_camera_model
struct PinholeCameraModel : public BaseCameraModel<PinholeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kPinhole, "PINHOLE", 2, 2, 0)
};

// Simple camera model with one focal length and one radial distortion
// parameter.
//
// This model is similar to the camera model that VisualSfM uses with the
// difference that the distortion here is applied to the projections and
// not to the measurements.
//
// Parameter list is expected in the following order:
//
//    f, cx, cy, k
//
struct SimpleRadialCameraModel
    : public BaseCameraModel<SimpleRadialCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kSimpleRadial, "SIMPLE_RADIAL", 1, 2, 1)
};

// Simple camera model with one focal length and two radial distortion
// parameters.
//
// This model is equivalent to the camera model that Bundler uses
// (except for an inverse z-axis in the camera coordinate system).
//
// Parameter list is expected in the following order:
//
//    f, cx, cy, k1, k2
//
struct RadialCameraModel : public BaseCameraModel<RadialCameraModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kRadial, "RADIAL", 1, 2, 2)
};

// OpenCV camera model.
//
// Based on the pinhole camera model. Additionally models radial and
// tangential distortion (up to 2nd degree of coefficients). Not suitable for
// large radial distortions of fish-eye cameras.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, k1, k2, p1, p2
//
// See
// http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html
struct OpenCVCameraModel : public BaseCameraModel<OpenCVCameraModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kOpenCV, "OPENCV", 2, 2, 4)
};

// OpenCV fish-eye camera model.
//
// Based on the pinhole camera model. Additionally models radial distortion
// (up to 4th degree of coefficients). Suitable for
// large radial distortions of fish-eye cameras.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, k1, k2, k3, k4
//
// See
// http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html
struct OpenCVFisheyeCameraModel
    : public BaseFisheyeCameraModel<OpenCVFisheyeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kOpenCVFisheye, "OPENCV_FISHEYE", 2, 2, 4)
  FISHEYE_CAMERA_MODEL_DEFINITIONS
};

// Full OpenCV camera model.
//
// Based on the pinhole camera model. Additionally models radial and
// tangential Distortion.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, k1, k2, p1, p2, k3, k4, k5, k6
//
// See
// http://docs.opencv.org/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html
struct FullOpenCVCameraModel : public BaseCameraModel<FullOpenCVCameraModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kFullOpenCV, "FULL_OPENCV", 2, 2, 8)
};

// FOV camera model.
//
// Based on the pinhole camera model. Additionally models radial distortion.
// This model is for example used by Project Tango for its equidistant
// calibration type.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, omega
//
// See:
// Frederic Devernay, Olivier Faugeras. Straight lines have to be straight:
// Automatic calibration and removal of distortion from scenes of structured
// environments. Machine vision and applications, 2001.
struct FOVCameraModel : public BaseCameraModel<FOVCameraModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kFOV, "FOV", 2, 2, 1)

  template <typename T>
  static void Undistortion(const T* extra_params, T u, T v, T* du, T* dv);
};

// Simple camera model with one focal length and one radial distortion
// parameter, suitable for fish-eye cameras.
//
// This model is equivalent to the OpenCVFisheyeCameraModel but has only one
// radial distortion coefficient.
//
// Parameter list is expected in the following order:
//
//    f, cx, cy, k
//
struct SimpleRadialFisheyeCameraModel
    : public BaseFisheyeCameraModel<SimpleRadialFisheyeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kSimpleRadialFisheye, "SIMPLE_RADIAL_FISHEYE", 1, 2, 1)
  FISHEYE_CAMERA_MODEL_DEFINITIONS
};

// Simple camera model with one focal length and two radial distortion
// parameters, suitable for fish-eye cameras.
//
// This model is equivalent to the OpenCVFisheyeCameraModel but has only two
// radial distortion coefficients.
//
// Parameter list is expected in the following order:
//
//    f, cx, cy, k1, k2
//
struct RadialFisheyeCameraModel
    : public BaseFisheyeCameraModel<RadialFisheyeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kRadialFisheye, "RADIAL_FISHEYE", 1, 2, 2)
  FISHEYE_CAMERA_MODEL_DEFINITIONS
};

// Camera model with radial and tangential distortion coefficients and
// additional coefficients accounting for thin-prism distortion.
//
// This camera model is described in
//
//    "Camera Calibration with Distortion Models and Accuracy Evaluation",
//    J Weng et al., TPAMI, 1992.
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, k1, k2, p1, p2, k3, k4, sx1, sy1
//
struct ThinPrismFisheyeCameraModel
    : public BaseFisheyeCameraModel<ThinPrismFisheyeCameraModel> {
  CAMERA_MODEL_DEFINITIONS(
      CameraModelId::kThinPrismFisheye, "THIN_PRISM_FISHEYE", 2, 2, 8)
  FISHEYE_CAMERA_MODEL_DEFINITIONS
};

// RadTanThinPrismFisheye Camera Model
//
// Camera model with radial and tangential distortion coefficients and
// additional coefficients accounting for thin-prism distortion.
//
// See
// https://facebookresearch.github.io/projectaria_tools/docs/tech_insights/camera_intrinsic_models#the-fisheyeradtanthinprism-fisheye624-model
//
// Parameter list is expected in the following order:
//
//    fx, fy, cx, cy, k0, k1, k2, k3, k4, k5, p0, p1, s0, s1, s2, s3
//
struct RadTanThinPrismFisheyeModel
    : public BaseFisheyeCameraModel<RadTanThinPrismFisheyeModel> {
  CAMERA_MODEL_DEFINITIONS(CameraModelId::kRadTanThinPrismFisheye,
                           "RAD_TAN_THIN_PRISM_FISHEYE",
                           2,
                           2,
                           12)
  FISHEYE_CAMERA_MODEL_DEFINITIONS
};

// Check whether camera model with given name or identifier exists.
bool ExistsCameraModelWithName(const std::string& model_name);
bool ExistsCameraModelWithId(CameraModelId model_id);

// Convert camera name to unique camera model identifier.
//
// @param name         Unique name of camera model.
//
// @return             Unique identifier of camera model.
CameraModelId CameraModelNameToId(const std::string& model_name);

// Convert camera model identifier to unique camera model name.
//
// @param model_id     Unique identifier of camera model.
//
// @return             Unique name of camera model.
const std::string& CameraModelIdToName(CameraModelId model_id);

// Initialize camera parameters using given image properties.
//
// Initializes all focal length parameters to the same given focal length and
// sets the principal point to the image center.
//
// @param model_id      Unique identifier of camera model.
// @param focal_length  Focal length, equal for all focal length parameters.
// @param width         Sensor width of the camera.
// @param height        Sensor height of the camera.
std::vector<double> CameraModelInitializeParams(CameraModelId model_id,
                                                double focal_length,
                                                size_t width,
                                                size_t height);

// Get human-readable information about the parameter vector order.
//
// @param model_id     Unique identifier of camera model.
const std::string& CameraModelParamsInfo(CameraModelId model_id);

// Get the indices of the parameter groups in the parameter vector.
//
// @param model_id     Unique identifier of camera model.
span<const size_t> CameraModelFocalLengthIdxs(CameraModelId model_id);
span<const size_t> CameraModelPrincipalPointIdxs(CameraModelId model_id);
span<const size_t> CameraModelExtraParamsIdxs(CameraModelId model_id);

// Get the total number of parameters of a camera model.
size_t CameraModelNumParams(CameraModelId model_id);

// Check whether parameters are valid, i.e. the parameter vector has
// the correct dimensions that match the specified camera model.
//
// @param model_id      Unique identifier of camera model.
// @param params        Array of camera parameters.
bool CameraModelVerifyParams(CameraModelId model_id,
                             const std::vector<double>& params);

// Check whether camera has bogus parameters.
//
// @param model_id                Unique identifier of camera model.
// @param params                  Array of camera parameters.
// @param width                   Sensor width of the camera.
// @param height                  Sensor height of the camera.
// @param min_focal_length_ratio  Minimum ratio of focal length over
//                                maximum sensor dimension.
// @param min_focal_length_ratio  Maximum ratio of focal length over
//                                maximum sensor dimension.
// @param max_extra_param         Maximum magnitude of each extra parameter.
bool CameraModelHasBogusParams(CameraModelId model_id,
                               const std::vector<double>& params,
                               size_t width,
                               size_t height,
                               double min_focal_length_ratio,
                               double max_focal_length_ratio,
                               double max_extra_param);

// Transform camera to image coordinates.
//
// This is the inverse of `CameraModelCamFromImg`.
//
// @param model_id     Unique model_id of camera model as defined in
//                     `CAMERA_MODEL_NAME_TO_CODE`.
// @param params       Array of camera parameters.
// @param u, v         Coordinates in camera system as (u, v, 1).
// @param x, y         Output image coordinates in pixels.
inline std::optional<Eigen::Vector2d> CameraModelImgFromCam(
    CameraModelId model_id,
    const std::vector<double>& params,
    const Eigen::Vector3d& uvw);

// Transform image to camera coordinates.
//
// This is the inverse of `CameraModelImgFromCam`.
//
// @param model_id      Unique identifier of camera model.
// @param params        Array of camera parameters.
// @param xy            Image coordinates in pixels.
//
// @return              Output ray in camera frame (u, v, w).
inline std::optional<Eigen::Vector2d> CameraModelCamFromImg(
    CameraModelId model_id,
    const std::vector<double>& params,
    const Eigen::Vector2d& xy);

// Convert pixel threshold in image plane to camera space by dividing
// the threshold through the mean focal length.
//
// @param model_id      Unique identifier of camera model.
// @param params        Array of camera parameters.
// @param threshold     Image space threshold in pixels.
//
// @return              Camera space threshold.
inline double CameraModelCamFromImgThreshold(CameraModelId model_id,
                                             const std::vector<double>& params,
                                             double threshold);

// Test if a camera model represents a fisheye camera.
//
// @param model_id      Unique identifier of camera model.
//
// @return              Whether it is a fisheye camera model.
inline bool CameraModelIsFisheye(CameraModelId model_id);

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// BaseCameraModel

template <typename CameraModel>
template <typename T>
bool BaseCameraModel<CameraModel>::HasBogusParams(
    const std::vector<T>& params,
    const size_t width,
    const size_t height,
    const T min_focal_length_ratio,
    const T max_focal_length_ratio,
    const T max_extra_param) {
  return HasBogusPrincipalPoint(params, width, height) ||
         HasBogusFocalLength(params,
                             width,
                             height,
                             min_focal_length_ratio,
                             max_focal_length_ratio) ||
         HasBogusExtraParams(params, max_extra_param);
}

template <typename CameraModel>
template <typename T>
bool BaseCameraModel<CameraModel>::HasBogusFocalLength(
    const std::vector<T>& params,
    const size_t width,
    const size_t height,
    const T min_focal_length_ratio,
    const T max_focal_length_ratio) {
  const T inv_max_size = 1.0 / std::max(width, height);
  for (const size_t idx : CameraModel::focal_length_idxs) {
    const T focal_length_ratio = params[idx] * inv_max_size;
    if (focal_length_ratio < min_focal_length_ratio ||
        focal_length_ratio > max_focal_length_ratio) {
      return true;
    }
  }

  return false;
}

template <typename CameraModel>
template <typename T>
bool BaseCameraModel<CameraModel>::HasBogusPrincipalPoint(
    const std::vector<T>& params, const size_t width, const size_t height) {
  const T cx = params[CameraModel::principal_point_idxs[0]];
  const T cy = params[CameraModel::principal_point_idxs[1]];
  return cx < 0 || cx > width || cy < 0 || cy > height;
}

template <typename CameraModel>
template <typename T>
bool BaseCameraModel<CameraModel>::HasBogusExtraParams(
    const std::vector<T>& params, const T max_extra_param) {
  for (const size_t idx : CameraModel::extra_params_idxs) {
    if (std::abs(params[idx]) > max_extra_param) {
      return true;
    }
  }

  return false;
}

template <typename CameraModel>
template <typename T>
T BaseCameraModel<CameraModel>::CamFromImgThreshold(const T* params,
                                                    const T threshold) {
  T mean_focal_length = 0;
  for (const size_t idx : CameraModel::focal_length_idxs) {
    mean_focal_length += params[idx];
  }
  mean_focal_length /= CameraModel::focal_length_idxs.size();
  return threshold / mean_focal_length;
}

template <typename CameraModel>
bool BaseCameraModel<CameraModel>::IterativeUndistortion(const double* params,
                                                         double* u,
                                                         double* v) {
  // Parameters for Newton iteration. 100 iterations should be enough for
  // complex camera models with higher order terms.
  constexpr size_t kNumIterations = 100;
  constexpr double kMinStepSquaredNorm = 1e-10;
  // Trust region: step_x.norm() <= max(x.norm() * kRelStepRadius, kStepRadius)
  constexpr double kRelStepRadius = 0.1;
  constexpr double kStepRadius = 0.1;

  Eigen::Matrix2d J;
  const Eigen::Vector2d x0(*u, *v);
  Eigen::Vector2d x(*u, *v);
  Eigen::Vector2d dx;

  ceres::Jet<double, 2> params_jet[CameraModel::num_extra_params];
  for (size_t i = 0; i < CameraModel::num_extra_params; ++i) {
    params_jet[i] = ceres::Jet<double, 2>(params[i]);
  }
  for (size_t i = 0; i < kNumIterations; ++i) {
    // Get Jacobian
    ceres::Jet<double, 2> x_jet[2];
    x_jet[0] = ceres::Jet<double, 2>(x(0), 0);
    x_jet[1] = ceres::Jet<double, 2>(x(1), 1);
    ceres::Jet<double, 2> dx_jet[2];
    CameraModel::Distortion(
        params_jet, x_jet[0], x_jet[1], &dx_jet[0], &dx_jet[1]);
    dx[0] = dx_jet[0].a;
    dx[1] = dx_jet[1].a;
    J(0, 0) = dx_jet[0].v[0] + 1;
    J(0, 1) = dx_jet[0].v[1];
    J(1, 0) = dx_jet[1].v[0];
    J(1, 1) = dx_jet[1].v[1] + 1;

    // Update
    Eigen::Vector2d step_x = J.partialPivLu().solve(x + dx - x0);
    const double radius_sqr =
        std::max(x.squaredNorm() * kRelStepRadius * kRelStepRadius,
                 kStepRadius * kStepRadius);
    const double step_x_norm_sqr = step_x.squaredNorm();
    if (step_x_norm_sqr > radius_sqr) {
      step_x *= std::sqrt(radius_sqr / step_x_norm_sqr);
    }
    x -= step_x;
    if (step_x.squaredNorm() < kMinStepSquaredNorm) {
      *u = x(0);
      *v = x(1);
      return true;
    }
  }

  *u = x(0);
  *v = x(1);

  return false;
}

////////////////////////////////////////////////////////////////////////////////
// SimplePinholeCameraModel

std::string SimplePinholeCameraModel::InitializeParamsInfo() {
  return "f, cx, cy";
}

std::array<size_t, 1> SimplePinholeCameraModel::InitializeFocalLengthIdxs() {
  return {0};
}

std::array<size_t, 2> SimplePinholeCameraModel::InitializePrincipalPointIdxs() {
  return {1, 2};
}

std::array<size_t, 0> SimplePinholeCameraModel::InitializeExtraParamsIdxs() {
  return {};
}

std::vector<double> SimplePinholeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, width / 2.0, height / 2.0};
}

template <typename T>
bool SimplePinholeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  // No Distortion

  // Transform to image coordinates
  *x = f * u / w + c1;
  *y = f * v / w + c2;

  return true;
}

bool SimplePinholeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f = params[0];
  const double c1 = params[1];
  const double c2 = params[2];

  *u = (x - c1) / f;
  *v = (y - c2) / f;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// PinholeCameraModel

std::string PinholeCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy";
}

std::array<size_t, 2> PinholeCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2> PinholeCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 0> PinholeCameraModel::InitializeExtraParamsIdxs() {
  return {};
}

std::vector<double> PinholeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, focal_length, width / 2.0, height / 2.0};
}

template <typename T>
bool PinholeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  // No Distortion

  // Transform to image coordinates
  *x = f1 * u / w + c1;
  *y = f2 * v / w + c2;

  return true;
}

bool PinholeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f1 = params[0];
  const double f2 = params[1];
  const double c1 = params[2];
  const double c2 = params[3];

  *u = (x - c1) / f1;
  *v = (y - c2) / f2;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// SimpleRadialCameraModel

std::string SimpleRadialCameraModel::InitializeParamsInfo() {
  return "f, cx, cy, k";
}

std::array<size_t, 1> SimpleRadialCameraModel::InitializeFocalLengthIdxs() {
  return {0};
}

std::array<size_t, 2> SimpleRadialCameraModel::InitializePrincipalPointIdxs() {
  return {1, 2};
}

std::array<size_t, 1> SimpleRadialCameraModel::InitializeExtraParamsIdxs() {
  return {3};
}

std::vector<double> SimpleRadialCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, width / 2.0, height / 2.0, 0};
}

template <typename T>
bool SimpleRadialCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  u /= w;
  v /= w;

  // Distortion
  T du, dv;
  Distortion(&params[3], u, v, &du, &dv);
  *x = u + du;
  *y = v + dv;

  // Transform to image coordinates
  *x = f * *x + c1;
  *y = f * *y + c2;

  return true;
}

bool SimpleRadialCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f = params[0];
  const double c1 = params[1];
  const double c2 = params[2];

  // Lift points to normalized plane
  *u = (x - c1) / f;
  *v = (y - c2) / f;

  return IterativeUndistortion(&params[3], u, v);
}

template <typename T>
void SimpleRadialCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k = extra_params[0];

  const T u2 = u * u;
  const T v2 = v * v;
  const T r2 = u2 + v2;
  const T radial = k * r2;
  *du = u * radial;
  *dv = v * radial;
}

////////////////////////////////////////////////////////////////////////////////
// RadialCameraModel

std::string RadialCameraModel::InitializeParamsInfo() {
  return "f, cx, cy, k1, k2";
}

std::array<size_t, 1> RadialCameraModel::InitializeFocalLengthIdxs() {
  return {0};
}

std::array<size_t, 2> RadialCameraModel::InitializePrincipalPointIdxs() {
  return {1, 2};
}

std::array<size_t, 2> RadialCameraModel::InitializeExtraParamsIdxs() {
  return {3, 4};
}

std::vector<double> RadialCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, width / 2.0, height / 2.0, 0, 0};
}

template <typename T>
bool RadialCameraModel::ImgFromCam(const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  u /= w;
  v /= w;

  // Distortion
  T du, dv;
  Distortion(&params[3], u, v, &du, &dv);
  *x = u + du;
  *y = v + dv;

  // Transform to image coordinates
  *x = f * *x + c1;
  *y = f * *y + c2;

  return true;
}

bool RadialCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f = params[0];
  const double c1 = params[1];
  const double c2 = params[2];

  // Lift points to normalized plane
  *u = (x - c1) / f;
  *v = (y - c2) / f;

  return IterativeUndistortion(&params[3], u, v);
}

template <typename T>
void RadialCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];

  const T u2 = u * u;
  const T v2 = v * v;
  const T r2 = u2 + v2;
  const T radial = k1 * r2 + k2 * r2 * r2;
  *du = u * radial;
  *dv = v * radial;
}

////////////////////////////////////////////////////////////////////////////////
// OpenCVCameraModel

std::string OpenCVCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, k1, k2, p1, p2";
}

std::array<size_t, 2> OpenCVCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2> OpenCVCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 4> OpenCVCameraModel::InitializeExtraParamsIdxs() {
  return {4, 5, 6, 7};
}

std::vector<double> OpenCVCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, focal_length, width / 2.0, height / 2.0, 0, 0, 0, 0};
}

template <typename T>
bool OpenCVCameraModel::ImgFromCam(const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  u /= w;
  v /= w;

  // Distortion
  T du, dv;
  Distortion(&params[4], u, v, &du, &dv);
  *x = u + du;
  *y = v + dv;

  // Transform to image coordinates
  *x = f1 * *x + c1;
  *y = f2 * *y + c2;

  return true;
}

bool OpenCVCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f1 = params[0];
  const double f2 = params[1];
  const double c1 = params[2];
  const double c2 = params[3];

  // Lift points to normalized plane
  *u = (x - c1) / f1;
  *v = (y - c2) / f2;

  return IterativeUndistortion(&params[4], u, v);
}

template <typename T>
void OpenCVCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];
  const T p1 = extra_params[2];
  const T p2 = extra_params[3];

  const T u2 = u * u;
  const T uv = u * v;
  const T v2 = v * v;
  const T r2 = u2 + v2;
  const T radial = k1 * r2 + k2 * r2 * r2;
  *du = u * radial + T(2) * p1 * uv + p2 * (r2 + T(2) * u2);
  *dv = v * radial + T(2) * p2 * uv + p1 * (r2 + T(2) * v2);
}

////////////////////////////////////////////////////////////////////////////////
// OpenCVFisheyeCameraModel

std::string OpenCVFisheyeCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, k1, k2, k3, k4";
}

std::array<size_t, 2> OpenCVFisheyeCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2> OpenCVFisheyeCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 4> OpenCVFisheyeCameraModel::InitializeExtraParamsIdxs() {
  return {4, 5, 6, 7};
}

std::vector<double> OpenCVFisheyeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, focal_length, width / 2.0, height / 2.0, 0, 0, 0, 0};
}

template <typename T>
void OpenCVFisheyeCameraModel::ImgFromFisheye(
    const T* params, T uu, T vv, T* x, T* y) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *x = f1 * uu + c1;
  *y = f2 * vv + c2;
}

template <typename T>
void OpenCVFisheyeCameraModel::FisheyeFromImg(
    const T* params, T x, T y, T* uu, T* vv) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *uu = (x - c1) / f1;
  *vv = (y - c2) / f2;
}

template <typename T>
bool OpenCVFisheyeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  u /= w;
  v /= w;
  T uu, vv;
  FisheyeFromNormal(u, v, &uu, &vv);

  // Distortion
  T duu, dvv;
  Distortion(&params[4], uu, vv, &duu, &dvv);

  // Transform to image coordinates
  ImgFromFisheye(params, uu + duu, vv + dvv, x, y);

  return true;
}

bool OpenCVFisheyeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  double uu, vv;
  FisheyeFromImg(params, x, y, &uu, &vv);
  if (!IterativeUndistortion(&params[4], &uu, &vv)) {
    return false;
  }
  NormalFromFisheye(uu, vv, u, v);
  return true;
}

template <typename T>
void OpenCVFisheyeCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];
  const T k3 = extra_params[2];
  const T k4 = extra_params[3];

  const T theta2 = u * u + v * v;
  const T theta4 = theta2 * theta2;
  const T theta6 = theta4 * theta2;
  const T theta8 = theta4 * theta4;
  const T radial = k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8;
  *du = u * radial;
  *dv = v * radial;
}

////////////////////////////////////////////////////////////////////////////////
// FullOpenCVCameraModel

std::string FullOpenCVCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, k1, k2, p1, p2, k3, k4, k5, k6";
}

std::array<size_t, 2> FullOpenCVCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2> FullOpenCVCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 8> FullOpenCVCameraModel::InitializeExtraParamsIdxs() {
  return {4, 5, 6, 7, 8, 9, 10, 11};
}

std::vector<double> FullOpenCVCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length,
          focal_length,
          width / 2.0,
          height / 2.0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0};
}

template <typename T>
bool FullOpenCVCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  u /= w;
  v /= w;

  // Distortion
  T du, dv;
  Distortion(&params[4], u, v, &du, &dv);
  *x = u + du;
  *y = v + dv;

  // Transform to image coordinates
  *x = f1 * *x + c1;
  *y = f2 * *y + c2;

  return true;
}

bool FullOpenCVCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f1 = params[0];
  const double f2 = params[1];
  const double c1 = params[2];
  const double c2 = params[3];

  // Lift points to normalized plane
  *u = (x - c1) / f1;
  *v = (y - c2) / f2;

  return IterativeUndistortion(&params[4], u, v);
}

template <typename T>
void FullOpenCVCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];
  const T p1 = extra_params[2];
  const T p2 = extra_params[3];
  const T k3 = extra_params[4];
  const T k4 = extra_params[5];
  const T k5 = extra_params[6];
  const T k6 = extra_params[7];

  const T u2 = u * u;
  const T uv = u * v;
  const T v2 = v * v;
  const T r2 = u2 + v2;
  const T r4 = r2 * r2;
  const T r6 = r4 * r2;
  const T radial = (T(1) + k1 * r2 + k2 * r4 + k3 * r6) /
                   (T(1) + k4 * r2 + k5 * r4 + k6 * r6);
  *du = u * radial + T(2) * p1 * uv + p2 * (r2 + T(2) * u2) - u;
  *dv = v * radial + T(2) * p2 * uv + p1 * (r2 + T(2) * v2) - v;
}

////////////////////////////////////////////////////////////////////////////////
// FOVCameraModel

std::string FOVCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, omega";
}

std::array<size_t, 2> FOVCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2> FOVCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 1> FOVCameraModel::InitializeExtraParamsIdxs() {
  return {4};
}

std::vector<double> FOVCameraModel::InitializeParams(const double focal_length,
                                                     const size_t width,
                                                     const size_t height) {
  return {focal_length, focal_length, width / 2.0, height / 2.0, 1e-2};
}

template <typename T>
bool FOVCameraModel::ImgFromCam(const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  u /= w;
  v /= w;

  // Distortion
  Distortion(&params[4], u, v, x, y);

  // Transform to image coordinates
  *x = f1 * *x + c1;
  *y = f2 * *y + c2;

  return true;
}

bool FOVCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  const double f1 = params[0];
  const double f2 = params[1];
  const double c1 = params[2];
  const double c2 = params[3];

  // Lift points to normalized plane
  const double uu = (x - c1) / f1;
  const double vv = (y - c2) / f2;

  // Undistortion
  Undistortion(&params[4], uu, vv, u, v);

  return true;
}

template <typename T>
void FOVCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T omega = extra_params[0];

  // Chosen arbitrarily.
  const T kEpsilon = T(1e-4);

  const T radius2 = u * u + v * v;
  const T omega2 = omega * omega;

  T factor;
  if (omega2 < kEpsilon) {
    // Derivation of this case with Matlab:
    // syms radius omega;
    // factor(radius) = atan(radius * 2 * tan(omega / 2)) / ...
    //                  (radius * omega);
    // simplify(taylor(factor, omega, 'order', 3))
    factor = (omega2 * radius2) / T(3) - omega2 / T(12) + T(1);
  } else if (radius2 < kEpsilon) {
    // Derivation of this case with Matlab:
    // syms radius omega;
    // factor(radius) = atan(radius * 2 * tan(omega / 2)) / ...
    //                  (radius * omega);
    // simplify(taylor(factor, radius, 'order', 3))
    const T tan_half_omega = ceres::tan(omega / T(2));
    factor = (T(-2) * tan_half_omega *
              (T(4) * radius2 * tan_half_omega * tan_half_omega - T(3))) /
             (T(3) * omega);
  } else {
    const T radius = ceres::sqrt(radius2);
    const T numerator = ceres::atan(radius * T(2) * ceres::tan(omega / T(2)));
    factor = numerator / (radius * omega);
  }

  *du = u * factor;
  *dv = v * factor;
}

template <typename T>
void FOVCameraModel::Undistortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  T omega = extra_params[0];

  // Chosen arbitrarily.
  const T kEpsilon = T(1e-4);

  const T radius2 = u * u + v * v;
  const T omega2 = omega * omega;

  T factor;
  if (omega2 < kEpsilon) {
    // Derivation of this case with Matlab:
    // syms radius omega;
    // factor(radius) = tan(radius * omega) / ...
    //                  (radius * 2*tan(omega/2));
    // simplify(taylor(factor, omega, 'order', 3))
    factor = (omega2 * radius2) / T(3) - omega2 / T(12) + T(1);
  } else if (radius2 < kEpsilon) {
    // Derivation of this case with Matlab:
    // syms radius omega;
    // factor(radius) = tan(radius * omega) / ...
    //                  (radius * 2*tan(omega/2));
    // simplify(taylor(factor, radius, 'order', 3))
    factor = (omega * (omega * omega * radius2 + T(3))) /
             (T(6) * ceres::tan(omega / T(2)));
  } else {
    const T radius = ceres::sqrt(radius2);
    const T numerator = ceres::tan(radius * omega);
    factor = numerator / (radius * T(2) * ceres::tan(omega / T(2)));
  }

  *du = u * factor;
  *dv = v * factor;
}

////////////////////////////////////////////////////////////////////////////////
// SimpleRadialFisheyeCameraModel

std::string SimpleRadialFisheyeCameraModel::InitializeParamsInfo() {
  return "f, cx, cy, k";
}

std::array<size_t, 1>
SimpleRadialFisheyeCameraModel::InitializeFocalLengthIdxs() {
  return {0};
}

std::array<size_t, 2>
SimpleRadialFisheyeCameraModel::InitializePrincipalPointIdxs() {
  return {1, 2};
}

std::array<size_t, 1>
SimpleRadialFisheyeCameraModel::InitializeExtraParamsIdxs() {
  return {3};
}

std::vector<double> SimpleRadialFisheyeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, width / 2.0, height / 2.0, 0};
}

template <typename T>
void SimpleRadialFisheyeCameraModel::ImgFromFisheye(
    const T* params, T uu, T vv, T* x, T* y) {
  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  *x = f * uu + c1;
  *y = f * vv + c2;
}

template <typename T>
void SimpleRadialFisheyeCameraModel::FisheyeFromImg(
    const T* params, T x, T y, T* uu, T* vv) {
  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  *uu = (x - c1) / f;
  *vv = (y - c2) / f;
}

template <typename T>
bool SimpleRadialFisheyeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  u /= w;
  v /= w;
  T uu, vv;
  FisheyeFromNormal(u, v, &uu, &vv);

  // Distortion
  T duu, dvv;
  Distortion(&params[3], uu, vv, &duu, &dvv);

  // Transform to image coordinates
  ImgFromFisheye(params, uu + duu, vv + dvv, x, y);

  return true;
}

bool SimpleRadialFisheyeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  double uu, vv;
  FisheyeFromImg(params, x, y, &uu, &vv);
  if (!IterativeUndistortion(&params[3], &uu, &vv)) {
    return false;
  }
  NormalFromFisheye(uu, vv, u, v);
  return true;
}

template <typename T>
void SimpleRadialFisheyeCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k = extra_params[0];

  const T theta2 = u * u + v * v;
  const T radial = k * theta2;
  *du = u * radial;
  *dv = v * radial;
}

////////////////////////////////////////////////////////////////////////////////
// RadialFisheyeCameraModel

std::string RadialFisheyeCameraModel::InitializeParamsInfo() {
  return "f, cx, cy, k1, k2";
}

std::array<size_t, 1> RadialFisheyeCameraModel::InitializeFocalLengthIdxs() {
  return {0};
}

std::array<size_t, 2> RadialFisheyeCameraModel::InitializePrincipalPointIdxs() {
  return {1, 2};
}

std::array<size_t, 2> RadialFisheyeCameraModel::InitializeExtraParamsIdxs() {
  return {3, 4};
}

std::vector<double> RadialFisheyeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length, width / 2.0, height / 2.0, 0, 0};
}

template <typename T>
void RadialFisheyeCameraModel::ImgFromFisheye(
    const T* params, T uu, T vv, T* x, T* y) {
  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  *x = f * uu + c1;
  *y = f * vv + c2;
}

template <typename T>
void RadialFisheyeCameraModel::FisheyeFromImg(
    const T* params, T x, T y, T* uu, T* vv) {
  const T f = params[0];
  const T c1 = params[1];
  const T c2 = params[2];

  *uu = (x - c1) / f;
  *vv = (y - c2) / f;
}

template <typename T>
bool RadialFisheyeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  u /= w;
  v /= w;
  T uu, vv;
  FisheyeFromNormal(u, v, &uu, &vv);

  // Distortion
  T duu, dvv;
  Distortion(&params[3], uu, vv, &duu, &dvv);

  // Transform to image coordinates
  ImgFromFisheye(params, uu + duu, vv + dvv, x, y);

  return true;
}

bool RadialFisheyeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  double uu, vv;
  FisheyeFromImg(params, x, y, &uu, &vv);
  if (!IterativeUndistortion(&params[3], &uu, &vv)) {
    return false;
  }
  NormalFromFisheye(uu, vv, u, v);
  return true;
}

template <typename T>
void RadialFisheyeCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];

  const T theta2 = u * u + v * v;
  const T theta4 = theta2 * theta2;
  const T radial = k1 * theta2 + k2 * theta4;
  *du = u * radial;
  *dv = v * radial;
}

////////////////////////////////////////////////////////////////////////////////
// ThinPrismFisheyeCameraModel

std::string ThinPrismFisheyeCameraModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, k1, k2, p1, p2, k3, k4, sx1, sy1";
}

std::array<size_t, 2> ThinPrismFisheyeCameraModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2>
ThinPrismFisheyeCameraModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 8> ThinPrismFisheyeCameraModel::InitializeExtraParamsIdxs() {
  return {4, 5, 6, 7, 8, 9, 10, 11};
}

std::vector<double> ThinPrismFisheyeCameraModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  return {focal_length,
          focal_length,
          width / 2.0,
          height / 2.0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0};
}

template <typename T>
void ThinPrismFisheyeCameraModel::ImgFromFisheye(
    const T* params, T uu, T vv, T* x, T* y) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *x = f1 * uu + c1;
  *y = f2 * vv + c2;
}

template <typename T>
void ThinPrismFisheyeCameraModel::FisheyeFromImg(
    const T* params, T x, T y, T* uu, T* vv) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *uu = (x - c1) / f1;
  *vv = (y - c2) / f2;
}

template <typename T>
bool ThinPrismFisheyeCameraModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  u /= w;
  v /= w;
  T uu, vv;
  FisheyeFromNormal(u, v, &uu, &vv);

  // Distortion
  T duu, dvv;
  Distortion(&params[4], uu, vv, &duu, &dvv);

  // Transform to image coordinates
  ImgFromFisheye(params, uu + duu, vv + dvv, x, y);

  return true;
}

bool ThinPrismFisheyeCameraModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  double uu, vv;
  FisheyeFromImg(params, x, y, &uu, &vv);
  if (!IterativeUndistortion(&params[4], &uu, &vv)) {
    return false;
  }
  NormalFromFisheye(uu, vv, u, v);
  return true;
}

template <typename T>
void ThinPrismFisheyeCameraModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  const T k1 = extra_params[0];
  const T k2 = extra_params[1];
  const T p1 = extra_params[2];
  const T p2 = extra_params[3];
  const T k3 = extra_params[4];
  const T k4 = extra_params[5];
  const T sx1 = extra_params[6];
  const T sy1 = extra_params[7];

  const T u2 = u * u;
  const T uv = u * v;
  const T v2 = v * v;
  const T r2 = u2 + v2;
  const T r4 = r2 * r2;
  const T r6 = r4 * r2;
  const T r8 = r6 * r2;
  const T radial = k1 * r2 + k2 * r4 + k3 * r6 + k4 * r8;
  *du = u * radial + T(2) * p1 * uv + p2 * (r2 + T(2) * u2) + sx1 * r2;
  *dv = v * radial + T(2) * p2 * uv + p1 * (r2 + T(2) * v2) + sy1 * r2;
}

////////////////////////////////////////////////////////////////////////////////
// RadTanThinPrismFisheyeModel

std::string RadTanThinPrismFisheyeModel::InitializeParamsInfo() {
  return "fx, fy, cx, cy, k0, k1, k2, k3, k4, k5, p0, p1, s0, s1, s2, s3";
}

std::array<size_t, 2> RadTanThinPrismFisheyeModel::InitializeFocalLengthIdxs() {
  return {0, 1};
}

std::array<size_t, 2>
RadTanThinPrismFisheyeModel::InitializePrincipalPointIdxs() {
  return {2, 3};
}

std::array<size_t, 12>
RadTanThinPrismFisheyeModel::InitializeExtraParamsIdxs() {
  return {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
}

std::vector<double> RadTanThinPrismFisheyeModel::InitializeParams(
    const double focal_length, const size_t width, const size_t height) {
  std::vector<double> params(RadTanThinPrismFisheyeModel::num_params);
  params[0] = focal_length;
  params[1] = focal_length;
  params[2] = width / 2.0;
  params[3] = height / 2.0;

  for (size_t i = 4; i < RadTanThinPrismFisheyeModel::num_params; ++i) {
    params[i] = 0;
  }

  return params;
}

template <typename T>
void RadTanThinPrismFisheyeModel::ImgFromFisheye(
    const T* params, T uu, T vv, T* x, T* y) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *x = f1 * uu + c1;
  *y = f2 * vv + c2;
}

template <typename T>
void RadTanThinPrismFisheyeModel::FisheyeFromImg(
    const T* params, T x, T y, T* uu, T* vv) {
  const T f1 = params[0];
  const T f2 = params[1];
  const T c1 = params[2];
  const T c2 = params[3];

  *uu = (x - c1) / f1;
  *vv = (y - c2) / f2;
}

template <typename T>
bool RadTanThinPrismFisheyeModel::ImgFromCam(
    const T* params, T u, T v, T w, T* x, T* y) {
  if (w < std::numeric_limits<T>::epsilon()) {
    return false;
  }

  u /= w;
  v /= w;
  T uu, vv;
  FisheyeFromNormal(u, v, &uu, &vv);

  T duu, dvv;
  Distortion(&params[4], uu, vv, &duu, &dvv);
  ImgFromFisheye(params, uu + duu, vv + dvv, x, y);

  return true;
}

bool RadTanThinPrismFisheyeModel::CamFromImg(
    const double* params, double x, double y, double* u, double* v) {
  double uu, vv;
  FisheyeFromImg(params, x, y, &uu, &vv);
  if (!IterativeUndistortion(&params[4], &uu, &vv)) {
    return false;
  }
  NormalFromFisheye(uu, vv, u, v);
  return true;
}

template <typename T>
void RadTanThinPrismFisheyeModel::Distortion(
    const T* extra_params, const T u, const T v, T* du, T* dv) {
  constexpr int kNumRadialParams = 6;
  const T radial_coeffs[kNumRadialParams] = {extra_params[0],
                                             extra_params[1],
                                             extra_params[2],
                                             extra_params[3],
                                             extra_params[4],
                                             extra_params[5]};

  const T p0 = extra_params[6];
  const T p1 = extra_params[7];
  const T s0 = extra_params[8];
  const T s1 = extra_params[9];
  const T s2 = extra_params[10];
  const T s3 = extra_params[11];

  const T theta2 = u * u + v * v;
  T th_radial = T(1);
  T theta_power = T(1);
  for (int i = 0; i < kNumRadialParams; ++i) {
    theta_power *= theta2;
    th_radial += radial_coeffs[i] * theta_power;
  }

  const T x = th_radial * u;
  const T y = th_radial * v;

  const T x2 = x * x;
  const T y2 = y * y;
  const T xy = x * y;
  const T r2 = x2 + y2;
  const T r4 = r2 * r2;

  const T dx_tang = T(2) * p1 * xy + p0 * (r2 + T(2) * x2);
  const T dy_tang = T(2) * p0 * xy + p1 * (r2 + T(2) * y2);

  const T dx_tp = s0 * r2 + s1 * r4;
  const T dy_tp = s2 * r2 + s3 * r4;

  const T x_distorted = x + dx_tang + dx_tp;
  const T y_distorted = y + dy_tang + dy_tp;

  *du = x_distorted - u;
  *dv = y_distorted - v;
}

std::optional<Eigen::Vector2d> CameraModelImgFromCam(
    const CameraModelId model_id,
    const std::vector<double>& params,
    const Eigen::Vector3d& uvw) {
  Eigen::Vector2d xy;
  switch (model_id) {
#define CAMERA_MODEL_CASE(CameraModel)                                     \
  case CameraModel::model_id:                                              \
    if (CameraModel::ImgFromCam(                                           \
            params.data(), uvw.x(), uvw.y(), uvw.z(), &xy.x(), &xy.y())) { \
      return xy;                                                           \
    }                                                                      \
    break;

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }
  return std::nullopt;
}

std::optional<Eigen::Vector2d> CameraModelCamFromImg(
    const CameraModelId model_id,
    const std::vector<double>& params,
    const Eigen::Vector2d& xy) {
  Eigen::Vector2d uv;
  switch (model_id) {
#define CAMERA_MODEL_CASE(CameraModel)                          \
  case CameraModel::model_id:                                   \
    if (CameraModel::CamFromImg(                                \
            params.data(), xy.x(), xy.y(), &uv.x(), &uv.y())) { \
      return uv;                                                \
    }                                                           \
    break;

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }
  return std::nullopt;
}

double CameraModelCamFromImgThreshold(const CameraModelId model_id,
                                      const std::vector<double>& params,
                                      const double threshold) {
  switch (model_id) {
#define CAMERA_MODEL_CASE(CameraModel)                                 \
  case CameraModel::model_id:                                          \
    return CameraModel::CamFromImgThreshold(params.data(), threshold); \
    break;

    CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
  }

  return -1;
}

bool CameraModelIsFisheye(const CameraModelId model_id) {
  switch (model_id) {
#define CAMERA_MODEL_CASE(CameraModel) case CameraModel::model_id:

    FISHEYE_CAMERA_MODEL_CASES
    return true;
    default:
      return false;

#undef CAMERA_MODEL_CASE
  }

  return false;
}

}  // namespace colmap
