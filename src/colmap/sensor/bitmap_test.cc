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

#include "colmap/sensor/bitmap.h"

#include "colmap/util/testing.h"

#include <FreeImage.h>
#include <gtest/gtest.h>

namespace colmap {
namespace {

TEST(Bitmap, BitmapColorEmpty) {
  BitmapColor<uint8_t> color;
  EXPECT_EQ(color.r, 0);
  EXPECT_EQ(color.g, 0);
  EXPECT_EQ(color.b, 0);
  EXPECT_EQ(color, BitmapColor<uint8_t>(0));
  EXPECT_EQ(color, BitmapColor<uint8_t>(0, 0, 0));
}

TEST(Bitmap, BitmapGrayColor) {
  BitmapColor<uint8_t> color(5);
  EXPECT_EQ(color.r, 5);
  EXPECT_EQ(color.g, 5);
  EXPECT_EQ(color.b, 5);
}

TEST(Bitmap, BitmapColorCast) {
  BitmapColor<float> color1(1.1f, 2.9f, -3.0f);
  BitmapColor<uint8_t> color2 = color1.Cast<uint8_t>();
  EXPECT_EQ(color2.r, 1);
  EXPECT_EQ(color2.g, 3);
  EXPECT_EQ(color2.b, 0);
}

TEST(Bitmap, Empty) {
  Bitmap bitmap;
  EXPECT_EQ(bitmap.Width(), 0);
  EXPECT_EQ(bitmap.Height(), 0);
  EXPECT_EQ(bitmap.Channels(), 0);
  EXPECT_FALSE(bitmap.IsRGB());
  EXPECT_FALSE(bitmap.IsGrey());
}

TEST(Bitmap, Print) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  std::ostringstream stream;
  stream << bitmap;
  EXPECT_EQ(stream.str(), "Bitmap(width=100, height=100, channels=3)");
}

TEST(Bitmap, AllocateRGB) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  EXPECT_EQ(bitmap.Width(), 100);
  EXPECT_EQ(bitmap.Height(), 100);
  EXPECT_EQ(bitmap.Channels(), 3);
  EXPECT_TRUE(bitmap.IsRGB());
  EXPECT_FALSE(bitmap.IsGrey());
}

TEST(Bitmap, AllocateGrey) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, false);
  EXPECT_EQ(bitmap.Width(), 100);
  EXPECT_EQ(bitmap.Height(), 100);
  EXPECT_EQ(bitmap.Channels(), 1);
  EXPECT_FALSE(bitmap.IsRGB());
  EXPECT_TRUE(bitmap.IsGrey());
}

TEST(Bitmap, Deallocate) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, false);
  bitmap.Deallocate();
  EXPECT_EQ(bitmap.Width(), 0);
  EXPECT_EQ(bitmap.Height(), 0);
  EXPECT_EQ(bitmap.Channels(), 0);
  EXPECT_EQ(bitmap.NumBytes(), 0);
  EXPECT_FALSE(bitmap.IsRGB());
  EXPECT_FALSE(bitmap.IsGrey());
}

TEST(Bitmap, MoveConstruct) {
  Bitmap bitmap;
  bitmap.Allocate(2, 1, true);
  const auto* data = bitmap.Data();
  Bitmap moved_bitmap(std::move(bitmap));
  EXPECT_EQ(moved_bitmap.Width(), 2);
  EXPECT_EQ(moved_bitmap.Height(), 1);
  EXPECT_EQ(moved_bitmap.Channels(), 3);
  EXPECT_EQ(moved_bitmap.Data(), data);
  // NOLINTBEGIN(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
  EXPECT_EQ(bitmap.Width(), 0);
  EXPECT_EQ(bitmap.Height(), 0);
  EXPECT_EQ(bitmap.Channels(), 0);
  EXPECT_EQ(bitmap.NumBytes(), 0);
  EXPECT_EQ(bitmap.Data(), nullptr);
  // NOLINTEND(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
}

TEST(Bitmap, MoveAssign) {
  Bitmap bitmap;
  bitmap.Allocate(2, 1, true);
  const auto* data = bitmap.Data();
  Bitmap moved_bitmap = std::move(bitmap);
  EXPECT_EQ(moved_bitmap.Width(), 2);
  EXPECT_EQ(moved_bitmap.Height(), 1);
  EXPECT_EQ(moved_bitmap.Channels(), 3);
  EXPECT_EQ(moved_bitmap.Data(), data);
  // NOLINTBEGIN(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
  EXPECT_EQ(bitmap.Width(), 0);
  EXPECT_EQ(bitmap.Height(), 0);
  EXPECT_EQ(bitmap.Channels(), 0);
  EXPECT_EQ(bitmap.NumBytes(), 0);
  EXPECT_EQ(bitmap.Data(), nullptr);
  // NOLINTEND(bugprone-use-after-move,clang-analyzer-cplusplus.Move)
}

TEST(Bitmap, BitsPerPixel) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  EXPECT_EQ(bitmap.BitsPerPixel(), 24);
  bitmap.Allocate(100, 100, false);
  EXPECT_EQ(bitmap.BitsPerPixel(), 8);
}

TEST(Bitmap, NumBytes) {
  Bitmap bitmap;
  EXPECT_EQ(bitmap.NumBytes(), 0);
  bitmap.Allocate(100, 100, true);
  EXPECT_EQ(bitmap.NumBytes(), 3 * 100 * 100);
  bitmap.Allocate(100, 100, false);
  EXPECT_EQ(bitmap.NumBytes(), 100 * 100);
}

TEST(Bitmap, ConvertToRowMajorArrayRGB) {
  Bitmap bitmap;
  bitmap.Allocate(2, 2, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  const std::vector<uint8_t> array = bitmap.ConvertToRowMajorArray();
  ASSERT_EQ(array.size(), 12);
  EXPECT_EQ(array[0], 0);
  EXPECT_EQ(array[1], 0);
  EXPECT_EQ(array[2], 0);
  EXPECT_EQ(array[3], 0);
  EXPECT_EQ(array[4], 0);
  EXPECT_EQ(array[5], 2);
  EXPECT_EQ(array[6], 0);
  EXPECT_EQ(array[7], 0);
  EXPECT_EQ(array[8], 1);
  EXPECT_EQ(array[9], 0);
  EXPECT_EQ(array[10], 0);
  EXPECT_EQ(array[11], 3);
}

TEST(Bitmap, ConvertToRowMajorArrayGrey) {
  Bitmap bitmap;
  bitmap.Allocate(2, 2, false);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  const std::vector<uint8_t> array = bitmap.ConvertToRowMajorArray();
  ASSERT_EQ(array.size(), 4);
  EXPECT_EQ(array[0], 0);
  EXPECT_EQ(array[1], 2);
  EXPECT_EQ(array[2], 1);
  EXPECT_EQ(array[3], 3);
}

TEST(Bitmap, ConvertToColMajorArrayRGB) {
  Bitmap bitmap;
  bitmap.Allocate(2, 2, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  const std::vector<uint8_t> array = bitmap.ConvertToColMajorArray();
  ASSERT_EQ(array.size(), 12);
  EXPECT_EQ(array[0], 0);
  EXPECT_EQ(array[1], 0);
  EXPECT_EQ(array[2], 0);
  EXPECT_EQ(array[3], 0);
  EXPECT_EQ(array[4], 0);
  EXPECT_EQ(array[5], 0);
  EXPECT_EQ(array[6], 0);
  EXPECT_EQ(array[7], 0);
  EXPECT_EQ(array[8], 0);
  EXPECT_EQ(array[9], 1);
  EXPECT_EQ(array[10], 2);
  EXPECT_EQ(array[11], 3);
}

TEST(Bitmap, ConvertToColMajorArrayGrey) {
  Bitmap bitmap;
  bitmap.Allocate(2, 2, false);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  const std::vector<uint8_t> array = bitmap.ConvertToColMajorArray();
  ASSERT_EQ(array.size(), 4);
  EXPECT_EQ(array[0], 0);
  EXPECT_EQ(array[1], 1);
  EXPECT_EQ(array[2], 2);
  EXPECT_EQ(array[3], 3);
}

TEST(Bitmap, ConvertToFromRawBitsGrey) {
  Bitmap bitmap;
  bitmap.Allocate(3, 2, false);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3));

  std::vector<uint8_t> raw_bits = bitmap.ConvertToRawBits();
  ASSERT_EQ(raw_bits.size(), bitmap.Pitch() * bitmap.Height());

  const std::vector<uint8_t> raw_bits_copy = raw_bits;
  Bitmap bitmap_copy = Bitmap::ConvertFromRawBits(raw_bits.data(),
                                                  bitmap.Pitch(),
                                                  bitmap.Width(),
                                                  bitmap.Height(),
                                                  /*rgb=*/false);
  EXPECT_EQ(bitmap.Width(), bitmap_copy.Width());
  EXPECT_EQ(bitmap.Height(), bitmap_copy.Height());
  EXPECT_EQ(bitmap.Channels(), bitmap_copy.Channels());
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(5));
  bitmap_copy.SetPixel(0, 1, BitmapColor<uint8_t>(5));
  EXPECT_EQ(raw_bits_copy, raw_bits);
  EXPECT_EQ(bitmap.ConvertToRowMajorArray(),
            bitmap_copy.ConvertToRowMajorArray());
}

TEST(Bitmap, ConvertToFromRawBitsRGB) {
  Bitmap bitmap;
  bitmap.Allocate(3, 2, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));

  std::vector<uint8_t> raw_bits = bitmap.ConvertToRawBits();
  ASSERT_EQ(raw_bits.size(), bitmap.Pitch() * bitmap.Height() * 3);

  const std::vector<uint8_t> raw_bits_copy = raw_bits;
  Bitmap bitmap_copy = Bitmap::ConvertFromRawBits(raw_bits.data(),
                                                  bitmap.Pitch(),
                                                  bitmap.Width(),
                                                  bitmap.Height(),
                                                  /*rgb=*/true);
  EXPECT_EQ(bitmap.Width(), bitmap_copy.Width());
  EXPECT_EQ(bitmap.Height(), bitmap_copy.Height());
  EXPECT_EQ(bitmap.Channels(), bitmap_copy.Channels());
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(5, 0, 0));
  bitmap_copy.SetPixel(0, 1, BitmapColor<uint8_t>(5, 0, 0));
  EXPECT_EQ(raw_bits_copy, raw_bits);
  EXPECT_EQ(bitmap.ConvertToRowMajorArray(),
            bitmap_copy.ConvertToRowMajorArray());
}

TEST(Bitmap, GetAndSetPixelRGB) {
  Bitmap bitmap;
  bitmap.Allocate(1, 1, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(1, 2, 3));
  BitmapColor<uint8_t> color;
  EXPECT_TRUE(bitmap.GetPixel(0, 0, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(1, 2, 3));
}

TEST(Bitmap, GetAndSetPixelGrey) {
  Bitmap bitmap;
  bitmap.Allocate(1, 1, false);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 2, 3));
  BitmapColor<uint8_t> color;
  EXPECT_TRUE(bitmap.GetPixel(0, 0, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(1, 2, 3));
  EXPECT_TRUE(bitmap.GetPixel(0, 0, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(1, 0, 0));
}

TEST(Bitmap, GetScanlineRGB) {
  Bitmap bitmap;
  bitmap.Allocate(3, 3, true);
  bitmap.Fill(BitmapColor<uint8_t>(1, 2, 3));
  for (size_t r = 0; r < 3; ++r) {
    const uint8_t* scanline = bitmap.GetScanline(r);
    for (size_t c = 0; c < 3; ++c) {
      BitmapColor<uint8_t> color;
      EXPECT_TRUE(bitmap.GetPixel(r, c, &color));
      EXPECT_EQ(scanline[c * 3 + FI_RGBA_RED], color.r);
      EXPECT_EQ(scanline[c * 3 + FI_RGBA_GREEN], color.g);
      EXPECT_EQ(scanline[c * 3 + FI_RGBA_BLUE], color.b);
    }
  }
}

TEST(Bitmap, GetScanlineGrey) {
  Bitmap bitmap;
  bitmap.Allocate(3, 3, false);
  bitmap.Fill(BitmapColor<uint8_t>(1, 2, 3));
  for (size_t r = 0; r < 3; ++r) {
    const uint8_t* scanline = bitmap.GetScanline(r);
    for (size_t c = 0; c < 3; ++c) {
      BitmapColor<uint8_t> color;
      EXPECT_TRUE(bitmap.GetPixel(r, c, &color));
      EXPECT_EQ(scanline[c], color.r);
    }
  }
}

TEST(Bitmap, Fill) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  bitmap.Fill(BitmapColor<uint8_t>(1, 2, 3));
  for (int y = 0; y < bitmap.Height(); ++y) {
    for (int x = 0; x < bitmap.Width(); ++x) {
      BitmapColor<uint8_t> color;
      EXPECT_TRUE(bitmap.GetPixel(x, y, &color));
      EXPECT_EQ(color, BitmapColor<uint8_t>(1, 2, 3));
    }
  }
}

TEST(Bitmap, InterpolateNearestNeighbor) {
  Bitmap bitmap;
  bitmap.Allocate(11, 11, true);
  bitmap.Fill(BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(5, 5, BitmapColor<uint8_t>(1, 2, 3));
  BitmapColor<uint8_t> color;
  EXPECT_TRUE(bitmap.InterpolateNearestNeighbor(5, 5, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(1, 2, 3));
  EXPECT_TRUE(bitmap.InterpolateNearestNeighbor(5.4999, 5.4999, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(1, 2, 3));
  EXPECT_TRUE(bitmap.InterpolateNearestNeighbor(5.5, 5.5, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(0, 0, 0));
  EXPECT_TRUE(bitmap.InterpolateNearestNeighbor(4.5, 5.4999, &color));
  EXPECT_EQ(color, BitmapColor<uint8_t>(1, 2, 3));
}

TEST(Bitmap, InterpolateBilinear) {
  Bitmap bitmap;
  bitmap.Allocate(11, 11, true);
  bitmap.Fill(BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(5, 5, BitmapColor<uint8_t>(1, 2, 3));
  BitmapColor<float> color;
  EXPECT_TRUE(bitmap.InterpolateBilinear(5, 5, &color));
  EXPECT_EQ(color, BitmapColor<float>(1, 2, 3));
  EXPECT_TRUE(bitmap.InterpolateBilinear(5.5, 5, &color));
  EXPECT_EQ(color, BitmapColor<float>(0.5, 1, 1.5));
  EXPECT_TRUE(bitmap.InterpolateBilinear(5.5, 5.5, &color));
  EXPECT_EQ(color, BitmapColor<float>(0.25, 0.5, 0.75));
}

TEST(Bitmap, SmoothRGB) {
  Bitmap bitmap;
  bitmap.Allocate(50, 50, true);
  for (int x = 0; x < 50; ++x) {
    for (int y = 0; y < 50; ++y) {
      bitmap.SetPixel(
          x, y, BitmapColor<uint8_t>(y * 50 + x, y * 50 + x, y * 50 + x));
    }
  }
  bitmap.Smooth(1, 1);
  EXPECT_EQ(bitmap.Width(), 50);
  EXPECT_EQ(bitmap.Height(), 50);
  EXPECT_EQ(bitmap.Channels(), 3);
  for (int x = 0; x < 50; ++x) {
    for (int y = 0; y < 50; ++y) {
      BitmapColor<uint8_t> color;
      EXPECT_TRUE(bitmap.GetPixel(x, y, &color));
      EXPECT_EQ(color.r, color.g);
      EXPECT_EQ(color.r, color.b);
    }
  }
}

TEST(Bitmap, SmoothGrey) {
  Bitmap bitmap;
  bitmap.Allocate(50, 50, false);
  for (int x = 0; x < 50; ++x) {
    for (int y = 0; y < 50; ++y) {
      bitmap.SetPixel(
          x, y, BitmapColor<uint8_t>(y * 50 + x, y * 50 + x, y * 50 + x));
    }
  }
  bitmap.Smooth(1, 1);
  EXPECT_EQ(bitmap.Width(), 50);
  EXPECT_EQ(bitmap.Height(), 50);
  EXPECT_EQ(bitmap.Channels(), 1);
}

TEST(Bitmap, RescaleRGB) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  Bitmap bitmap1 = bitmap.Clone();
  bitmap1.Rescale(50, 25);
  EXPECT_EQ(bitmap1.Width(), 50);
  EXPECT_EQ(bitmap1.Height(), 25);
  EXPECT_EQ(bitmap1.Channels(), 3);
  Bitmap bitmap2 = bitmap.Clone();
  bitmap2.Rescale(150, 20);
  EXPECT_EQ(bitmap2.Width(), 150);
  EXPECT_EQ(bitmap2.Height(), 20);
  EXPECT_EQ(bitmap2.Channels(), 3);
}

TEST(Bitmap, RescaleGrey) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, false);
  Bitmap bitmap1 = bitmap.Clone();
  bitmap1.Rescale(50, 25);
  EXPECT_EQ(bitmap1.Width(), 50);
  EXPECT_EQ(bitmap1.Height(), 25);
  EXPECT_EQ(bitmap1.Channels(), 1);
  Bitmap bitmap2 = bitmap.Clone();
  bitmap2.Rescale(150, 20);
  EXPECT_EQ(bitmap2.Width(), 150);
  EXPECT_EQ(bitmap2.Height(), 20);
  EXPECT_EQ(bitmap2.Channels(), 1);
}

TEST(Bitmap, Clone) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  const Bitmap cloned_bitmap = bitmap.Clone();
  EXPECT_EQ(cloned_bitmap.Width(), 100);
  EXPECT_EQ(cloned_bitmap.Height(), 100);
  EXPECT_EQ(cloned_bitmap.Channels(), 3);
  EXPECT_NE(bitmap.Data(), cloned_bitmap.Data());
}

TEST(Bitmap, CloneAsRGB) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, false);
  const Bitmap cloned_bitmap = bitmap.CloneAsRGB();
  EXPECT_EQ(cloned_bitmap.Width(), 100);
  EXPECT_EQ(cloned_bitmap.Height(), 100);
  EXPECT_EQ(cloned_bitmap.Channels(), 3);
  EXPECT_NE(bitmap.Data(), cloned_bitmap.Data());
}

TEST(Bitmap, CloneAsGrey) {
  Bitmap bitmap;
  bitmap.Allocate(100, 100, true);
  const Bitmap cloned_bitmap = bitmap.CloneAsGrey();
  EXPECT_EQ(cloned_bitmap.Width(), 100);
  EXPECT_EQ(cloned_bitmap.Height(), 100);
  EXPECT_EQ(cloned_bitmap.Channels(), 1);
  EXPECT_NE(bitmap.Data(), cloned_bitmap.Data());
}

TEST(Bitmap, ReadWriteAsRGB) {
  Bitmap bitmap;
  bitmap.Allocate(2, 3, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  bitmap.SetPixel(0, 2, BitmapColor<uint8_t>(4, 2, 0));
  bitmap.SetPixel(1, 2, BitmapColor<uint8_t>(5, 2, 1));

  const std::string test_dir = CreateTestDir();
  const std::string filename = test_dir + "/bitmap.png";

  EXPECT_TRUE(bitmap.Write(filename));

  Bitmap read_bitmap;

  // Allocate bitmap with different size to test read overwrites existing data.
  read_bitmap.Allocate(bitmap.Width() + 1, bitmap.Height() + 2, true);

  EXPECT_TRUE(read_bitmap.Read(filename));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 3);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 24);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.ConvertToRowMajorArray());

  EXPECT_TRUE(read_bitmap.Read(filename, /*as_rgb=*/false));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 1);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 8);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.CloneAsGrey().ConvertToRowMajorArray());
}

TEST(Bitmap, ReadWriteAsGrey) {
  Bitmap bitmap;
  bitmap.Allocate(2, 3, false);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3));
  bitmap.SetPixel(0, 2, BitmapColor<uint8_t>(4));
  bitmap.SetPixel(1, 2, BitmapColor<uint8_t>(5));

  const std::string test_dir = CreateTestDir();
  const std::string filename = test_dir + "/bitmap.png";

  EXPECT_TRUE(bitmap.Write(filename));

  Bitmap read_bitmap;

  // Allocate bitmap with different size to test read overwrites existing data.
  read_bitmap.Allocate(bitmap.Width() + 1, bitmap.Height() + 2, true);

  EXPECT_TRUE(read_bitmap.Read(filename));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 3);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 24);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.CloneAsRGB().ConvertToRowMajorArray());

  EXPECT_TRUE(read_bitmap.Read(filename, /*as_rgb=*/false));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 1);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 8);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.ConvertToRowMajorArray());
}

TEST(Bitmap, ReadRGB16AsGrey) {
  Bitmap bitmap;
  bitmap.Allocate(2, 3, true);
  bitmap.SetPixel(0, 0, BitmapColor<uint8_t>(0, 0, 0));
  bitmap.SetPixel(0, 1, BitmapColor<uint8_t>(1, 0, 0));
  bitmap.SetPixel(1, 0, BitmapColor<uint8_t>(2, 0, 0));
  bitmap.SetPixel(1, 1, BitmapColor<uint8_t>(3, 0, 0));
  bitmap.SetPixel(0, 2, BitmapColor<uint8_t>(4, 2, 0));
  bitmap.SetPixel(1, 2, BitmapColor<uint8_t>(5, 2, 1));

  const std::string test_dir = CreateTestDir();
  const std::string filename = test_dir + "/bitmap.png";

  // Bitmap class does not support 16 bit color depth
  FIBITMAP* converted_rgb16 = FreeImage_ConvertToType(bitmap.Data(), FIT_RGB16);
  EXPECT_TRUE(converted_rgb16);
  EXPECT_TRUE(FreeImage_Save(FIF_PNG, converted_rgb16, filename.c_str()));
  FreeImage_Unload(converted_rgb16);

  // Assert the file was written correctly with 16 bit color depth
  FIBITMAP* written_image = FreeImage_Load(FIF_PNG, filename.c_str());
  EXPECT_TRUE(written_image);
  EXPECT_EQ(FreeImage_GetBPP(written_image), 48);
  FreeImage_Unload(written_image);

  Bitmap read_bitmap;

  // Allocate bitmap with different size to test read overwrites existing data.
  read_bitmap.Allocate(bitmap.Width() + 1, bitmap.Height() + 2, true);

  EXPECT_TRUE(read_bitmap.Read(filename));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 3);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 24);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.ConvertToRowMajorArray());

  EXPECT_TRUE(read_bitmap.Read(filename, /*as_rgb=*/false));
  EXPECT_EQ(read_bitmap.Width(), bitmap.Width());
  EXPECT_EQ(read_bitmap.Height(), bitmap.Height());
  EXPECT_EQ(read_bitmap.Channels(), 1);
  EXPECT_EQ(read_bitmap.BitsPerPixel(), 8);
  EXPECT_EQ(read_bitmap.ConvertToRowMajorArray(),
            bitmap.CloneAsGrey().ConvertToRowMajorArray());
}

}  // namespace
}  // namespace colmap
