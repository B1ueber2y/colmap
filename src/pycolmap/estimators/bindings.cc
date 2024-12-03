#include <pybind11/pybind11.h>

namespace py = pybind11;

void BindAbsolutePoseEstimator(py::module& m);
void BindAlignmentEstimator(py::module& m);
void BindBundleAdjuster(py::module& m);
void BindCovarianceEstimator(py::module& m);
void BindEssentialMatrixEstimator(py::module& m);
void BindFundamentalMatrixEstimator(py::module& m);
void BindGeneralizedAbsolutePoseEstimator(py::module& m);
void BindHomographyMatrixEstimator(py::module& m);
void BindSimilarityTransformEstimator(py::module& m);
void BindTriangulationEstimator(py::module& m);
void BindTwoViewGeometryEstimator(py::module& m);
#ifdef PYCOLMAP_PYCERES_ENABLED
void BindCostFunctions(py::module& m);
void BindManifold(py::module& m);
#endif

void BindEstimators(py::module& m) {
  BindAbsolutePoseEstimator(m);
  BindAlignmentEstimator(m);
  BindBundleAdjuster(m);
  BindCovarianceEstimator(m);
  BindEssentialMatrixEstimator(m);
  BindFundamentalMatrixEstimator(m);
  BindGeneralizedAbsolutePoseEstimator(m);
  BindHomographyMatrixEstimator(m);
  BindSimilarityTransformEstimator(m);
  BindTriangulationEstimator(m);
  BindTwoViewGeometryEstimator(m);
#ifdef PYCOLMAP_PYCERES_ENABLED
  BindCostFunctions(m);
  BindManifold(m);
#endif
}
