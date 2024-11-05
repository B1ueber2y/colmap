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

#include <vector>
#include <map>
#include <unordered_map>

namespace colmap {

struct ControlPoint {
  ControlPoint() {}
  ControlPoint(int sequence_id, int cp_id, std::string name, std::pair<time_t, time_t> timestamps): sequence_id(sequence_id), cp_id(cp_id), name(name), timestamps(timestamps) {}

  int sequence_id;
  int cp_id;
  std::string name;
  std::pair<time_t, time_t> timestamps; // (low, high)
};

struct Segment {
  Segment() {}
  Segment(int sequence_id, int cp_id_low, int cp_id_hight): sequence_id(sequence_id), cp_id_low(cp_id_low), cp_id_high(cp_id_high) {}

  int sequence_id;
  int cp_id_low;
  int cp_id_high;
}

class ControlPointSequence {
public:
  enum NodeType {
    CP = 0,
    SEGMENT = 1
  };

  ControlPointSequence() {}
  explicit ControlPointSequence(const std::vector<ControlPoint>& control_points);
  int sequence_id;
  std::vector<ControlPoint> control_points;
  std::vector<Segment> segments;

  // TODO: import images
  std::unordered_map<image_t, std::vector<NodeType, int>>; // image_id -> Node 
};

class CrossSequenceMatching {
  CrossSequenceMatches(seq_id_1, seq_id_2, matches): sequence_id_1(seq_id_1), sequence_id_2(seq_id_2), matches(matches) {}
  int sequence_id_1;
  int sequence_id_2;
  std::vector<std::pair<image_t, image_t>> matches;
};

class ControlPointSegmentGraph {
public:
  // TODO: match cps by name when importing sequences
  void ImportSequence(const ControlPointSequence* sequence);
  void ImportCrossSequenceMatching(const CrossSequenceMatching& matches);

  void GetNeighboringRanges(ControlPoint* cp, std::map<int, std::pair<time_t, time_t>>* ranges, int depth=3);
  void GetNeighboringRanges(Segment* segment, std::map<int, std::pair<time_t, time_t>>* ranges, int depth=3);

private:
  void AddControlPoint(ControlPoint* cp);
  void AddSegment(ControlPoint* cp);

  void AddEdge(ControlPoint* cp, Segment* segment);
  void AddEdge(ControlPoint* cp1, ControlPoint* cp2);
  void AddEdge(Segment* segment1, Segment* segment2);

  void GetNeighboringSegments(ControlPoint* cp, std::vector<Segment*>* neighboring_segments);
  void GetNeighboringSegments(Segment* segment, std::vector<Segment*>* neighboring_segments);
  void GetNeighboringControlPoints(ControlPoint* cp, std::vector<ControlPoint*>* neighboring_cps);
  void GetNeighboringControlPoints(Segment* segment, std::vector<ControlPoint*>* neighboring_cps);

  std::unordered_map<ControlPoint*, std::vector<ControlPoint*>> m_cp_to_cp;
  std::unordered_map<ControlPoint*, std::vector<Segment*>> m_cp_to_segment;
  std::unordered_map<Segment*, std::vector<ControlPoint*>> m_segment_to_cp;
  std::unordered_map<Segment*, std::vector<Segment*>> m_segment_to_segment;

  std::unordered_map<int, ControlPointSequence*> sequences;
};

}  // namespace colmap
