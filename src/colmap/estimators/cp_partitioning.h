// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
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

#include "colmap/util/types.h"

#include <map>
#include <set>
#include <vector>

namespace colmap {

typedef double timestamp_t;

struct ControlPoint {
  ControlPoint() {}
  ControlPoint(int sequence_id,
               int cp_id,
               std::string name,
               std::pair<timestamp_t, timestamp_t> timestamps)
      : sequence_id(sequence_id),
        cp_id(cp_id),
        name(name),
        timestamps(timestamps) {}
  bool operator<(const ControlPoint& other) const {
    if (sequence_id < other.sequence_id)
      return true;
    else if (sequence_id == other.sequence_id && cp_id < other.cp_id)
      return true;
    return false;
  };

  int sequence_id;
  int cp_id;
  std::string name;
  std::pair<timestamp_t, timestamp_t> timestamps;  // (low, high)
};

struct Segment {
  Segment() {}
  Segment(int sequence_id, int segment_id, int cp_id_low, int cp_id_high)
      : sequence_id(sequence_id),
        segment_id(segment_id),
        cp_id_low(cp_id_low),
        cp_id_high(cp_id_high) {}
  bool operator<(const Segment& other) const {
    if (sequence_id < other.sequence_id)
      return true;
    else if (sequence_id == other.sequence_id && segment_id < other.segment_id)
      return true;
    return false;
  }

  int sequence_id;
  int segment_id;
  int cp_id_low;
  int cp_id_high;
};

class ControlPointSequence {
 public:
  enum NodeType { CP = 0, SEGMENT = 1 };

  ControlPointSequence() {}
  ControlPointSequence(const std::vector<ControlPoint>& control_points);
  bool operator<(const ControlPointSequence& other) const {
    return sequence_id < other.sequence_id;
  }
  int sequence_id;
  std::vector<ControlPoint> control_points;
  std::vector<Segment> segments;

  void ImportImages(const std::map<image_t, timestamp_t>& images);

 private:
  std::map<image_t, std::pair<NodeType, int>> images;  // image_id -> Node
};

class SequenceMatching {
 public:
  SequenceMatching() {}
  SequenceMatching(int seq_id_1,
                   int seq_id_2,
                   const std::vector<std::pair<image_t, image_t>>& matches)
      : sequence_id_1(seq_id_1), sequence_id_2(seq_id_2), matches(matches) {}
  int sequence_id_1;
  int sequence_id_2;
  std::vector<std::pair<image_t, image_t>> matches;
};

using ControlPointIndex = std::pair<int, int>;
using SegmentIndex = std::pair<int, int>;
class ControlPointSegmentGraph {
 public:
  ControlPointSegmentGraph() {}
  // TODO: match cps by name when importing sequences
  void ImportSequence(ControlPointSequence* sequence);
  void ImportSequenceMatching(const SequenceMatching& matches);

  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      const ControlPoint& base_cp, int maxDepth = 3) const;
  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      const Segment& base_segment, int maxDepth = 3) const;

 private:
  void AddControlPoint(const ControlPoint& cp);
  bool HasControlPoint(const ControlPointIndex& cp_index) const;
  bool HasControlPoint(const ControlPoint& cp) const;
  void AddSegment(const Segment& segment);
  bool HasSegment(const SegmentIndex& segment_index) const;
  bool HasSegment(const Segment& segment) const;

  ControlPointIndex GetControlPointIndex(const ControlPoint& cp) const;
  ControlPoint GetControlPoint(const ControlPointIndex& cp_index) const;
  SegmentIndex GetSegmentIndex(const Segment& segment) const;
  Segment GetSegment(const SegmentIndex& segment_index) const;

  void AddEdge(const ControlPoint& cp, const Segment& segment);
  void AddEdge(const ControlPoint& cp1, const ControlPoint& cp2);
  void AddEdge(const Segment& segment1, const Segment& segment2);

  void GetNeighboringSegmentsFromCP(
      const ControlPointIndex& cp_index,
      std::vector<SegmentIndex>* neighboring_segments) const;
  void GetNeighboringSegmentsFromSegment(
      const SegmentIndex& segment_index,
      std::vector<SegmentIndex>* neighboring_segments) const;
  void GetNeighboringControlPointsFromCP(
      const ControlPointIndex& cp_index,
      std::vector<ControlPointIndex>* neighboring_cps) const;
  void GetNeighboringControlPointsFromSegment(
      const SegmentIndex& segment_index,
      std::vector<ControlPointIndex>* neighboring_cps) const;

  std::map<ControlPointIndex, std::set<ControlPointIndex>> m_cp_to_cp;
  std::map<ControlPointIndex, std::set<SegmentIndex>> m_cp_to_segment;
  std::map<SegmentIndex, std::set<ControlPointIndex>> m_segment_to_cp;
  std::map<SegmentIndex, std::set<SegmentIndex>> m_segment_to_segment;

  std::map<int, ControlPointSequence*> sequences_;
};

}  // namespace colmap
