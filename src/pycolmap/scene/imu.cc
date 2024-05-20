#include "colmap/sensor/imu.h"

#include "colmap/scene/imu.h"

#include "pycolmap/helpers.h"
#include "pycolmap/pybind11_extension.h"
#include "pycolmap/utils.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace colmap;
using namespace pybind11::literals;
namespace py = pybind11;

void BindImu(py::module& m) {
  py::class_<ImuCalibration>(m, "ImuCalibration")
      .def(py::init<>())
      .def_readwrite("acc_noise_density", &ImuCalibration::acc_noise_density)
      .def_readwrite("gyro_noise_density", &ImuCalibration::gyro_noise_density)
      .def_readwrite("acc_bias_random_walk_sigma",
                     &ImuCalibration::acc_bias_random_walk_sigma)
      .def_readwrite("gyro_bias_random_walk_sigma",
                     &ImuCalibration::gyro_bias_random_walk_sigma)
      .def_readwrite("acc_saturation_max", &ImuCalibration::acc_saturation_max)
      .def_readwrite("gyro_saturation_max",
                     &ImuCalibration::gyro_saturation_max)
      .def_readwrite("gravity_magnitude", &ImuCalibration::gravity_magnitude)
      .def_readwrite("imu_rate", &ImuCalibration::imu_rate);

  py::class_<ImuMeasurement>(m, "ImuMeasurement")
      .def(py::init<const double,
                    const Eigen::Vector3d&,
                    const Eigen::Vector3d&>())
      .def_readwrite("timestamp", &ImuMeasurement::timestamp)
      .def_readwrite("linear_acceleration",
                     &ImuMeasurement::linear_acceleration)
      .def_readwrite("angular_velocity", &ImuMeasurement::angular_velocity)
      .def("__repr__", [](const ImuMeasurement& m) {
        std::stringstream ss;
        ss << "ImuMeasurement("
           << "t=" << m.timestamp << ", "
           << "acc=[" << m.linear_acceleration.format(vec_fmt) << "], "
           << "gyro=[" << m.angular_velocity.format(vec_fmt) << "])";
        return ss.str();
      });

  py::class_<ImuMeasurements>(m, "ImuMeasurements")
      .def(py::init<>())
      .def(py::init<const std::vector<ImuMeasurement>&>())
      .def(py::init<const ImuMeasurements&>())
      .def("insert",
           py::overload_cast<const ImuMeasurement&>(&ImuMeasurements::insert))
      .def("insert",
           py::overload_cast<const std::vector<ImuMeasurement>&>(
               &ImuMeasurements::insert))
      .def("insert",
           py::overload_cast<const ImuMeasurements&>(&ImuMeasurements::insert))
      .def("remove", &ImuMeasurements::remove)
      .def("front", &ImuMeasurements::front)
      .def("back", &ImuMeasurements::back)
      .def("empty", &ImuMeasurements::empty)
      .def("clear", &ImuMeasurements::clear)
      .def("__len__", &ImuMeasurements::size)
      .def("__getitem__",
           [](const ImuMeasurements& ms, size_t idx) { return ms[idx]; })
      .def("get_measurements_contain_edge",
           &ImuMeasurements::GetMeasurementsContainEdge)
      .def_property_readonly("data", &ImuMeasurements::Data)
      .def("__repr__", [](const ImuMeasurements& ms) {
        std::stringstream ss;
        ss << "ImuMeasurements("
           << "n=" << ms.size() << ", "
           << "tmin=" << ms.front().timestamp << ", "
           << "tmax=" << ms.back().timestamp << ")";
        return ss.str();
      });

  py::class_<Imu>(m, "Imu")
      .def(py::init<>())
      .def_readwrite("imu_id", &Imu::camera_id)
      .def_readwrite("camera_id", &Imu::camera_id)
      .def_readwrite("cam_to_imu", &Imu::cam_to_imu);

  py::class_<ImuState>(m, "ImuState")
      .def(py::init<>())
      .def("set_velocity", &ImuState::SetVelocity)
      .def("set_acc_bias", &ImuState::SetAccBias)
      .def("set_gyro_bias", &ImuState::SetGyroBias)
      .def_property_readonly("data", &ImuState::Data)
      .def_property_readonly("velocity", &ImuState::Velocity)
      .def_property_readonly("velocity_ptr", &ImuState::VelocityPtr)
      .def_property_readonly("acc_bias", &ImuState::AccBias)
      .def_property_readonly("acc_bias_ptr", &ImuState::AccBiasPtr)
      .def_property_readonly("gyro_bias", &ImuState::GyroBias)
      .def_property_readonly("gyro_bias_ptr", &ImuState::GyroBiasPtr);
}