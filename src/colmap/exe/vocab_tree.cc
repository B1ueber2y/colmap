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

#include "colmap/exe/vocab_tree.h"

#include "colmap/controllers/feature_matching.h"
#include "colmap/controllers/option_manager.h"
#include "colmap/feature/sift.h"
#include "colmap/feature/utils.h"
#include "colmap/optim/random_sampler.h"
#include "colmap/retrieval/resources.h"
#include "colmap/retrieval/visual_index.h"
#include "colmap/scene/database.h"
#include "colmap/util/file.h"
#include "colmap/util/opengl_utils.h"

#include <numeric>

namespace colmap {
namespace {

// Loads descriptors for training from the database. Loads all descriptors from
// the database if max_num_images < 0, otherwise the descriptors of a random
// subset of images are selected.
retrieval::VisualIndex::Descriptors LoadRandomDatabaseDescriptors(
    const std::string& database_path, const int max_num_images) {
  Database database(database_path);
  DatabaseTransaction database_transaction(&database);

  const std::vector<Image> images = database.ReadAllImages();

  retrieval::VisualIndex::Descriptors descriptors;

  std::vector<size_t> image_idxs;
  size_t num_descriptors = 0;
  if (max_num_images < 0) {
    // All images in the database.
    image_idxs.resize(images.size());
    std::iota(image_idxs.begin(), image_idxs.end(), 0);
    num_descriptors = database.NumDescriptors();
  } else {
    // Random subset of images in the database.
    THROW_CHECK_LE(max_num_images, images.size());
    RandomSampler random_sampler(max_num_images);
    random_sampler.Initialize(images.size());
    random_sampler.Sample(&image_idxs);
    for (const size_t image_idx : image_idxs) {
      const auto& image = images.at(image_idx);
      num_descriptors += database.NumDescriptorsForImage(image.ImageId());
    }
  }

  descriptors.resize(num_descriptors, 128);

  size_t descriptor_row = 0;
  for (const size_t image_idx : image_idxs) {
    const auto& image = images.at(image_idx);
    const FeatureDescriptors image_descriptors =
        database.ReadDescriptors(image.ImageId());
    descriptors.block(descriptor_row, 0, image_descriptors.rows(), 128) =
        image_descriptors.cast<float>();
    descriptor_row += image_descriptors.rows();
  }

  THROW_CHECK_EQ(descriptor_row, num_descriptors);

  return descriptors;
}

std::vector<Image> ReadVocabTreeRetrievalImageList(const std::string& path,
                                                   Database* database) {
  std::vector<Image> images;
  if (path.empty()) {
    images.reserve(database->NumImages());
    for (const auto& image : database->ReadAllImages()) {
      images.push_back(image);
    }
  } else {
    DatabaseTransaction database_transaction(database);

    const std::vector<std::string> image_names = ReadTextFileLines(path);
    images.reserve(image_names.size());
    for (const auto& image_name : image_names) {
      Image image = database->ReadImageWithName(image_name).value();
      THROW_CHECK_NE(image.ImageId(), kInvalidImageId);
      images.push_back(std::move(image));
    }
  }
  return images;
}

}  // namespace

int RunVocabTreeBuilder(int argc, char** argv) {
  std::string vocab_tree_path = kDefaultVocabTreeUri;
  retrieval::VisualIndex::BuildOptions build_options;
  int max_num_images = -1;

  OptionManager options;
  options.AddDatabaseOptions();
  options.AddRequiredOption("vocab_tree_path", &vocab_tree_path);
  options.AddDefaultOption("num_visual_words", &build_options.num_visual_words);
  options.AddDefaultOption("num_iterations", &build_options.num_iterations);
  options.AddDefaultOption("num_checks", &build_options.num_checks);
  options.AddDefaultOption("num_threads", &build_options.num_threads);
  options.AddDefaultOption("num_rounds", &build_options.num_rounds);
  options.AddDefaultOption("max_num_images", &max_num_images);
  options.Parse(argc, argv);

  LOG(INFO) << "Loading descriptors...";
  const auto descriptors =
      LoadRandomDatabaseDescriptors(*options.database_path, max_num_images);
  LOG(INFO) << "=> Loaded a total of " << descriptors.rows() << " descriptors";
  THROW_CHECK_GT(descriptors.size(), 0);

  auto visual_index = retrieval::VisualIndex::Create();

  LOG(INFO) << "Building index for visual words...";
  // NOLINTNEXTLINE(clang-analyzer-optin.cplusplus.VirtualCall)
  visual_index->Build(build_options, descriptors);
  LOG(INFO) << "=> Quantized descriptor space using "
            << visual_index->NumVisualWords() << " visual words";

  LOG(INFO) << "Saving index to file...";
  visual_index->Write(vocab_tree_path);

  return EXIT_SUCCESS;
}

int RunVocabTreeRetriever(int argc, char** argv) {
  std::string vocab_tree_path = kDefaultVocabTreeUri;
  std::string database_image_list_path;
  std::string query_image_list_path;
  std::string output_index_path;
  retrieval::VisualIndex::QueryOptions query_options;
  retrieval::VisualIndex::IndexOptions index_options;
  int max_num_features = -1;

  OptionManager options;
  options.AddDatabaseOptions();
  options.AddRequiredOption("vocab_tree_path", &vocab_tree_path);
  options.AddDefaultOption("database_image_list_path",
                           &database_image_list_path);
  options.AddDefaultOption("query_image_list_path", &query_image_list_path);
  options.AddDefaultOption("output_index_path", &output_index_path);
  options.AddDefaultOption("num_images", &query_options.max_num_images);
  options.AddDefaultOption("num_neighbors", &query_options.num_neighbors);
  options.AddDefaultOption("num_checks", &query_options.num_checks);
  options.AddDefaultOption("num_threads", &query_options.num_threads);
  options.AddDefaultOption("num_images_after_verification",
                           &query_options.num_images_after_verification);
  options.AddDefaultOption("max_num_features", &max_num_features);
  options.Parse(argc, argv);

  index_options.num_threads = query_options.num_threads;

  auto visual_index = retrieval::VisualIndex::Read(vocab_tree_path);

  Database database(*options.database_path);

  const auto database_images =
      ReadVocabTreeRetrievalImageList(database_image_list_path, &database);
  const auto query_images =
      (!query_image_list_path.empty() || output_index_path.empty())
          ? ReadVocabTreeRetrievalImageList(query_image_list_path, &database)
          : std::vector<Image>();

  //////////////////////////////////////////////////////////////////////////////
  // Perform image indexing
  //////////////////////////////////////////////////////////////////////////////

  for (size_t i = 0; i < database_images.size(); ++i) {
    Timer timer;
    timer.Start();

    LOG(INFO) << StringPrintf(
                     "Indexing image [%d/%d]", i + 1, database_images.size())
              << std::flush;

    if (visual_index->IsImageIndexed(database_images[i].ImageId())) {
      continue;
    }

    FeatureKeypoints keypoints =
        database.ReadKeypoints(database_images[i].ImageId());
    FeatureDescriptors descriptors =
        database.ReadDescriptors(database_images[i].ImageId());
    if (max_num_features > 0 && descriptors.rows() > max_num_features) {
      ExtractTopScaleFeatures(&keypoints, &descriptors, max_num_features);
    }

    visual_index->Add(index_options,
                      database_images[i].ImageId(),
                      keypoints,
                      descriptors.cast<float>());

    LOG(INFO) << StringPrintf(" in %.3fs", timer.ElapsedSeconds());
  }

  // Compute the TF-IDF weights, etc.
  visual_index->Prepare();

  // Optionally save the indexing data for the database images (as well as the
  // original vocabulary tree data) to speed up future indexing.
  if (!output_index_path.empty()) {
    visual_index->Write(output_index_path);
  }

  if (query_images.empty()) {
    return EXIT_SUCCESS;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Perform image queries
  //////////////////////////////////////////////////////////////////////////////

  std::unordered_map<image_t, const Image*> image_id_to_image;
  image_id_to_image.reserve(database_images.size());
  for (const auto& image : database_images) {
    image_id_to_image.emplace(image.ImageId(), &image);
  }

  for (size_t i = 0; i < query_images.size(); ++i) {
    Timer timer;
    timer.Start();

    LOG(INFO) << StringPrintf("Querying for image %s [%d/%d]",
                              query_images[i].Name().c_str(),
                              i + 1,
                              query_images.size())
              << std::flush;

    auto keypoints = database.ReadKeypoints(query_images[i].ImageId());
    auto descriptors = database.ReadDescriptors(query_images[i].ImageId());
    if (max_num_features > 0 && descriptors.rows() > max_num_features) {
      ExtractTopScaleFeatures(&keypoints, &descriptors, max_num_features);
    }

    std::vector<retrieval::ImageScore> image_scores;
    visual_index->Query(
        query_options, keypoints, descriptors.cast<float>(), &image_scores);

    LOG(INFO) << StringPrintf(" in %.3fs", timer.ElapsedSeconds());
    for (const auto& image_score : image_scores) {
      const auto& image = *image_id_to_image.at(image_score.image_id);
      LOG(INFO) << StringPrintf("  image_id=%d, image_name=%s, score=%f",
                                image_score.image_id,
                                image.Name().c_str(),
                                image_score.score);
    }
  }

  return EXIT_SUCCESS;
}

}  // namespace colmap
