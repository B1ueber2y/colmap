#include "colmap/scene/track.h"

#include "colmap/util/logging.h"
#include "colmap/util/misc.h"
#include "colmap/util/types.h"

#include "pycolmap/helpers.h"

#include <memory>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

using namespace colmap;
using namespace pybind11::literals;
namespace py = pybind11;

void BindTrack(py::module& m) {
  py::class_<TrackElement, std::shared_ptr<TrackElement>> PyTrackElement(
      m, "TrackElement");
  PyTrackElement.def(py::init<>())
      .def(py::init<image_t, point2D_t>(), "image_id"_a, "point2D_idx"_a)
      .def_readwrite("image_id", &TrackElement::image_id)
      .def_readwrite("point2D_idx", &TrackElement::point2D_idx);
  MakeDataclass(PyTrackElement);

  py::class_<Track, std::shared_ptr<Track>> PyTrack(m, "Track");
  PyTrack.def(py::init<>())
      .def(py::init([](const std::vector<TrackElement>& elements) {
             auto track = std::make_shared<Track>();
             track->AddElements(elements);
             return track;
           }),
           "elements"_a)
      .def_property("elements",
                    py::overload_cast<>(&Track::Elements),
                    &Track::SetElements)
      .def("length", &Track::Length, "Track Length.")
      .def("add_element",
           py::overload_cast<image_t, point2D_t>(&Track::AddElement),
           "image_id"_a,
           "point2D_idx"_a,
           "Add an observation to the track.")
      .def("add_element",
           py::overload_cast<const TrackElement&>(&Track::AddElement),
           "element"_a)
      .def("add_elements",
           &Track::AddElements,
           "elements"_a,
           "Add multiple elements.")
      .def("delete_element",
           py::overload_cast<image_t, point2D_t>(&Track::DeleteElement),
           "image_id"_a,
           "point2D_idx"_a,
           "Delete observation from track.")
      .def("delete_element",
           py::overload_cast<size_t>(&Track::DeleteElement),
           "index"_a,
           "Remove TrackElement at index.");
  MakeDataclass(PyTrack);
}
