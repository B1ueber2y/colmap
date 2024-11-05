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
#include <unordered_map>
#include <vector>

namespace colmap {

typedef int32_t time_t;

struct ControlPoint {
  ControlPoint() {}
  ControlPoint(int sequence_id,
               int cp_id,
               std::string name,
               std::pair<time_t, time_t> timestamps)
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
  std::pair<time_t, time_t> timestamps;  // (low, high)
};

struct Segment {
  Segment() {}
  Segment(int sequence_id, int cp_id_low, int cp_id_high)
      : sequence_id(sequence_id),
        cp_id_low(cp_id_low),
        cp_id_high(cp_id_high) {}
  bool operator<(const Segment& other) const {
    if (sequence_id < other.sequence_id)
      return true;
    else if (sequence_id == other.sequence_id && cp_id_low < other.cp_id_low)
      return true;
    else if (sequence_id == other.sequence_id && cp_id_low == other.cp_id_low &&
             cp_id_high < other.cp_id_high)
      return true;
    return false;
  }

  int sequence_id;
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

  // TODO: import images
  std::unordered_map<image_t, std::pair<NodeType, int>>
      images;  // image_id -> Node
};

class CrossSequenceMatching {
  CrossSequenceMatching(int seq_id_1,
                        int seq_id_2,
                        const std::vector<std::pair<image_t, image_t>>& matches)
      : sequence_id_1(seq_id_1), sequence_id_2(seq_id_2), matches(matches) {}
  int sequence_id_1;
  int sequence_id_2;
  std::vector<std::pair<image_t, image_t>> matches;
};

class ControlPointSegmentGraph {
 public:
  ControlPointSegmentGraph() {}
  // TODO: match cps by name when importing sequences
  void ImportSequence(ControlPointSequence* sequence);
  void ImportCrossSequenceMatching(const CrossSequenceMatching& matches);

  std::map<int, std::pair<time_t, time_t>> GetNeighboringRanges(
      ControlPoint base_cp, int maxDepth = 3) const;
  std::map<int, std::pair<time_t, time_t>> GetNeighboringRanges(
      Segment base_segment, int maxDepth = 3) const;

 private:
  void AddControlPoint(ControlPoint cp);
  bool HasControlPoint(ControlPoint cp) const;
  void AddSegment(Segment segment);
  bool HasSegment(Segment segment) const;

  void AddEdge(ControlPoint* cp, Segment* segment);
  void AddEdge(ControlPoint* cp1, ControlPoint* cp2);
  void AddEdge(Segment* segment1, Segment* segment2);

  void GetNeighboringSegments(
      ControlPoint cp, std::vector<Segment*>* neighboring_segments) const;
  void GetNeighboringSegments(
      Segment segment, std::vector<Segment*>* neighboring_segments) const;
  void GetNeighboringControlPoints(
      ControlPoint cp, std::vector<ControlPoint*>* neighboring_cps) const;
  void GetNeighboringControlPoints(
      Segment segment, std::vector<ControlPoint*>* neighboring_cps) const;

  std::map<ControlPoint, std::vector<ControlPoint*>> m_cp_to_cp;
  std::map<ControlPoint, std::vector<Segment*>> m_cp_to_segment;
  std::map<Segment, std::vector<ControlPoint*>> m_segment_to_cp;
  std::map<Segment, std::vector<Segment*>> m_segment_to_segment;

  std::map<int, ControlPointSequence*> sequences_;
};

}  // namespace colmap
