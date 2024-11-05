#include "colmap/estimators/cp_partitioning.h"

#include "pycolmap/helpers.h"
#include "pycolmap/pybind11_extension.h"
#include "pycolmap/utils.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace colmap;
using namespace pybind11::literals;
namespace py = pybind11;

void BindCPPartitioning(py::module& m) {
  py::class_<ControlPoint> PyControlPoint(m, "ControlPoint");
  PyControlPoint.def(py::init<>())
    .def(py::init<int, int, std::string, std::pair<timestamp_t, timestamp_t>>())
    .def_readwrite("sequence_id", &ControlPoint::sequence_id)
    .def_readwrite("cp_id", &ControlPoint::cp_id)
    .def_readwrite("name", &ControlPoint::name)
    .def_readwrite("timestamps", &ControlPoint::timestamps);

  py::class_<Segment> PySegment(m, "Segment");
  PySegment.def(py::init<>())
    .def_readwrite("sequence_id", &Segment::sequence_id)
    .def_readwrite("cp_id_low", &Segment::cp_id_low)
    .def_readwrite("cp_id_high", &Segment::cp_id_high);

  py::class_<ControlPointSequence> PyControlPointSequence(m, "ControlPointSequence");
  PyControlPointSequence.def(py::init<>())
    .def(py::init<const std::vector<ControlPoint>&>())
    .def_readonly("sequence_id", &ControlPointSequence::sequence_id)
    .def_readonly("control_points", &ControlPointSequence::control_points)
    .def_readonly("segments", &ControlPointSequence::segments);

  py::class_<ControlPointSegmentGraph>(m, "ControlPointSegmentGraph")
    .def(py::init<>())
    .def("import_sequence", &ControlPointSegmentGraph::ImportSequence)
    .def("get_neighboring_ranges", py::overload_cast<ControlPoint, int>(&ControlPointSegmentGraph::GetNeighboringRanges, py::const_),
        py::arg("base_cp"), py::arg("max_depth") = 3)
    .def("get_neighboring_ranges", py::overload_cast<Segment, int>(&ControlPointSegmentGraph::GetNeighboringRanges, py::const_),
        py::arg("base_segment"), py::arg("max_depth") = 3);

}
