#include "colmap/estimators/bundle_adjustment.h"
#include "colmap/scene/reconstruction.h"
#include "colmap/scene/synthetic.h"

#include <benchmark/benchmark.h>

using namespace colmap;

class BM_BundleAdjustment : public benchmark::Fixture {
 public:
   void SetUp(::benchmark::State& state) {
     int input_size = state.range(0);
     int ratio = input_size / default_size;

     SyntheticDatasetOptions synthetic_dataset_options;
     synthetic_dataset_options.num_cameras = 3 * ratio;
     synthetic_dataset_options.num_images = 10 * ratio;
     synthetic_dataset_options.num_points3D = 200 * ratio;
     synthetic_dataset_options.point2D_stddev = 0.01;
     SynthesizeDataset(synthetic_dataset_options, reconstruction.get());
   
     BundleAdjustmentConfig config;
     for (const auto& [image_id, image] : reconstruction->Images()) {
       config.AddImage(image_id);
     }

     // Fix the Gauge by always setting at least 3 points as constant.
     CHECK_GT(reconstruction->NumPoints3D(), 3);
     int num_constant_points = 0;
     for (const auto& [point3D_id, _] : reconstruction->Points3D()) {
       if (++num_constant_points <= 3) {
         config.AddConstantPoint(point3D_id);
       }
     }
     bundle_adjuster = CreateDefaultBundleAdjuster(BundleAdjustmentOptions(), std::move(config), *reconstruction);
   }

   void TearDown(::benchmark::State& state) {
     reconstruction->TearDown();
   }

   std::shared_ptr<Reconstruction> reconstruction;
   const int default_size = 10;
   std::unique_ptr<BundleAdjuster> bundle_adjuster;
};

BENCHMARK_DEFINE_F(BM_BundleAdjustment, Run)(benchmark::State& state) {
  for (auto _: state) {
    auto summary = bundle_adjuster->Solve();
  }
}

BENCHMARK_REGISTER_F(BM_BundleAdjustment, Run)->Arg(10)->Arg(50);

BENCHMARK_MAIN();
