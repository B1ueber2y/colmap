using time_t = int64_t
using InertialBias = Eigen::Vector6d;

struct Gauge {
  std::optional<double> scale;
  std::optional<Eigen::Vector3d> gravity;
  std::optional<double> heading;
  std::optional<Eigen::Vector3d> origin;
};

using rig_t = uint64_t;
using frame_t = uint64_t;

enum class SensorType {
  Camera = 0,
  IMU = 1,
  Location = 2, // include GNSS, radios, compass, etc.
};

struct sensor_t {
  rig_t id;
  SensorType type;
};

struct Image {
  data_t data_id;
  time_t timestamp = -1;
  std::shared_ptr<CameraCalibration> calibration;
  std::shared_ptr<Frame> frame = nullptr; // nullptr if not registered

  bool IsRegistered() { return frame != nullptr; }

  data_t ImageId() { return data_id; }

};

struct Frame {
  frame_t frame_id;
  std::vector<data_t> data_id;
  rig_t3d rig_from_world;
  std::shared_ptr<RigCalibration> rig_calibration; // nullptr if no rig
  rig_t3d SensorFromWorld() {
    CHECK_EQ(rig_calibration, nullptr);
    return rig_from_world;
  }
  rig_t3d SensorFromWorld(sensor_t sensor_id) {
    if (rig_calibration != nullptr) {
      return rig_calibration.SensorFromRig(sensor_id) * rig_from_world;
    }
    else
      return rig_from_world;
  }
};

struct InertialFrame: public Frame {
  Eigen::Vector3d velocity;
};

struct RigCalibration {
  rig_t rig_id;
  std::unordered_map<sensor_t, rig_t3d> sensor_from_rig;
};

struct ImageGroup {
  time_t timestamp;
  std::shared_ptr<RigCalibration> rig_calibration; // nullptr if a single image
};
