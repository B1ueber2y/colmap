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
  using namespace colmap::cp_partitioning;
  py::class_<ControlPoint> PyControlPoint(m, "ControlPoint");
  PyControlPoint.def(py::init<>())
      .def(py::init<int,
                    int,
                    std::string,
                    std::pair<timestamp_t, timestamp_t>>(),
           "sequence_id"_a,
           "cp_id"_a,
           "name"_a,
           "timestamps"_a)
      .def_readwrite("sequence_id", &ControlPoint::sequence_id)
      .def_readwrite("cp_id", &ControlPoint::cp_id)
      .def_readwrite("name", &ControlPoint::name)
      .def_readwrite("timestamps", &ControlPoint::timestamps);

  py::class_<Segment> PySegment(m, "Segment");
  PySegment.def(py::init<>())
      .def_readwrite("sequence_id", &Segment::sequence_id)
      .def_readwrite("segment_id", &Segment::segment_id)
      .def_readwrite("cp_id_left", &Segment::cp_id_left)
      .def_readwrite("cp_id_right", &Segment::cp_id_right);

  py::class_<ControlPointSequence> PyControlPointSequence(
      m, "ControlPointSequence");
  PyControlPointSequence.def(py::init<>())
      .def(py::init<const std::vector<ControlPoint>&,
                    const std::pair<timestamp_t, timestamp_t>&>())
      .def("import_image_timestamps",
           &ControlPointSequence::ImportImageTimestamps)
      .def("get_index", &ControlPointSequence::GetIndex)
      .def("get_image_ids_inside_time_range",
           &ControlPointSequence::GetImageIdsInsideTimeRange)
      .def("get_image_ids_from_cp", &ControlPointSequence::GetImageIdsFromCP)
      .def("get_image_ids_from_segment",
           &ControlPointSequence::GetImageIdsFromSegment)
      .def("get_image_ids_from_elements",
           [](ControlPointSequence& self, const py::list& items) {
             std::vector<std::pair<NodeType, int>> nodes;
             for (auto& item : items) {
               if (py::isinstance<ControlPoint>(item)) {
                 nodes.push_back(std::make_pair(
                     NodeType::CP, item.cast<ControlPoint>().cp_id));
               } else if (py::isinstance<Segment>(item)) {
                 nodes.push_back(std::make_pair(
                     NodeType::SEGMENT, item.cast<Segment>().segment_id));
               } else {
                 throw py::type_error("Error! Unsupported type");
               }
             }
             return self.GetImageIdsFromNodeCollection(nodes);
           })
      .def_readonly("sequence_id", &ControlPointSequence::sequence_id)
      .def_readonly("control_points", &ControlPointSequence::control_points)
      .def_readonly("segments", &ControlPointSequence::segments)
      .def_readonly("time_ranges", &ControlPointSequence::time_ranges);

  py::class_<SequenceMatching> PySequenceMatching(m, "SequenceMatching");
  PySequenceMatching.def(py::init<>())
      .def(
          py::init<int, int, const std::vector<std::pair<image_t, image_t>>&>(),
          "sequence_id_1"_a,
          "sequence_id_2"_a,
          "matches"_a)
      .def_readonly("sequence_id_1", &SequenceMatching::sequence_id_1)
      .def_readonly("sequence_id_2", &SequenceMatching::sequence_id_2)
      .def_readonly("matches", &SequenceMatching::matches);

  py::class_<ControlPointSegmentGraph>(m, "ControlPointSegmentGraph")
      .def(py::init<>())
      .def_readonly("sequences", &ControlPointSegmentGraph::sequences)
      .def("num_control_points", &ControlPointSegmentGraph::NumControlPoints)
      .def("num_segments", &ControlPointSegmentGraph::NumSegments)
      .def("num_nodes", &ControlPointSegmentGraph::NumNodes)
      .def("num_edges", &ControlPointSegmentGraph::NumEdges)
      .def("import_sequence", &ControlPointSegmentGraph::ImportSequence)
      .def("import_sequence_matching",
           &ControlPointSegmentGraph::ImportSequenceMatching)
      .def("import_matching_from_reconstruction",
           &ControlPointSegmentGraph::ImportMatchingFromReconstruction,
           py::arg("reconstruction"),
           py::arg("min_num_shared_points") = 30)
      .def("get_neighboring_ranges",
           py::overload_cast<const ControlPoint&, int>(
               &ControlPointSegmentGraph::GetNeighboringRanges, py::const_),
           py::arg("base_cp"),
           py::arg("max_depth") = 3)
      .def("get_neighboring_ranges",
           py::overload_cast<const Segment&, int>(
               &ControlPointSegmentGraph::GetNeighboringRanges, py::const_),
           py::arg("base_segment"),
           py::arg("max_depth") = 3)
      .def(
          "get_neighboring_ranges",
          [](ControlPointSegmentGraph& self,
             const py::list& items,
             int maxDepth) {
            std::vector<Node> nodes;
            for (auto& item : items) {
              if (py::isinstance<ControlPoint>(item)) {
                nodes.push_back(self.GetNode(item.cast<ControlPoint>()));
              } else if (py::isinstance<Segment>(item)) {
                nodes.push_back(self.GetNode(item.cast<Segment>()));
              } else {
                throw py::type_error("Error! Unsupported type");
              }
            }
            return self.GetNeighboringRanges(nodes, maxDepth);
          },
          py::arg("base_elements"),
          py::arg("max_depth") = 3)
      .def("get_neighboring_image_ids",
           py::overload_cast<const ControlPoint&, int>(
               &ControlPointSegmentGraph::GetNeighboringImageIds, py::const_),
           py::arg("base_cp"),
           py::arg("max_depth") = 3)
      .def("get_neighboring_image_ids",
           py::overload_cast<const Segment&, int>(
               &ControlPointSegmentGraph::GetNeighboringImageIds, py::const_),
           py::arg("base_segment"),
           py::arg("max_depth") = 3)
      .def(
          "get_neighboring_image_ids",
          [](ControlPointSegmentGraph& self,
             const py::list& items,
             int maxDepth) {
            std::vector<Node> nodes;
            for (auto& item : items) {
              if (py::isinstance<ControlPoint>(item)) {
                nodes.push_back(self.GetNode(item.cast<ControlPoint>()));
              } else if (py::isinstance<Segment>(item)) {
                nodes.push_back(self.GetNode(item.cast<Segment>()));
              } else {
                throw py::type_error("Error! Unsupported type");
              }
            }
            return self.GetNeighboringImageIds(nodes, maxDepth);
          },
          py::arg("base_elements"),
          py::arg("max_depth") = 3);
}
