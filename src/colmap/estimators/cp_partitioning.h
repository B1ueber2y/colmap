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

#include "colmap/scene/reconstruction.h"
#include "colmap/util/types.h"

#include <map>
#include <queue>
#include <set>
#include <vector>

namespace colmap {

namespace cp_partitioning {

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

  int sequence_id;
  int cp_id;
  std::string name;
  std::pair<timestamp_t, timestamp_t> timestamps;  // (low, high)
};

struct Segment {
  Segment() {}
  Segment(int sequence_id, int segment_id, int cp_id_left, int cp_id_right)
      : sequence_id(sequence_id),
        segment_id(segment_id),
        cp_id_left(cp_id_left),
        cp_id_right(cp_id_right) {}

  int sequence_id;
  int segment_id;
  int cp_id_left;   // -1 means the segment is at thes tart of the session
  int cp_id_right;  // -1 means the segment is at the end of the session
};

enum NodeType { CP = 0, SEGMENT = 1 };
class ControlPointSequence {
 public:
  ControlPointSequence() {}
  ControlPointSequence(const std::vector<ControlPoint>& control_points,
                       const std::pair<timestamp_t, timestamp_t>& time_ranges);
  bool operator<(const ControlPointSequence& other) const {
    return sequence_id < other.sequence_id;
  }
  int sequence_id;
  std::vector<ControlPoint> control_points;
  std::vector<Segment> segments;
  std::pair<timestamp_t, timestamp_t> time_ranges;

  timestamp_t GetSegmentStartTime(const int segment_id) const;
  timestamp_t GetSegmentEndTime(const int segment_id) const;

  void ImportImageTimestamps(
      const std::map<image_t, timestamp_t>& images_timestamps);
  std::pair<NodeType, int> GetIndex(const image_t image_id) const;

  std::vector<image_t> GetImageIds() const;
  std::vector<image_t> GetImageIdsInsideTimeRange(
      const std::pair<timestamp_t, timestamp_t>& time_range) const;
  std::vector<image_t> GetImageIdsFromCP(const ControlPoint& cp) const;
  std::vector<image_t> GetImageIdsFromSegment(const Segment& segment) const;
  std::vector<image_t> GetImageIdsFromNodeCollection(
      const std::vector<std::pair<NodeType, int>>& indexes) const;

 private:
  std::map<image_t, timestamp_t> image_timestamps_;
  std::map<image_t, std::pair<NodeType, int>>
      images_;  // image_id -> cp / segment
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

using NodeIndex = std::pair<int, int>;
using Node = std::pair<NodeType, NodeIndex>;
class ControlPointSegmentGraph {
 public:
  ControlPointSegmentGraph() {}
  // interfaces
  void ImportSequence(const ControlPointSequence& sequence);
  void ImportSequenceMatching(const SequenceMatching& matches);
  void ImportMatchingFromReconstruction(const Reconstruction& reconstruction,
                                        int min_num_shared_point = 30);

  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      const ControlPoint& base_cp, int maxDepth = 3) const;
  std::vector<image_t> GetNeighboringImageIds(const ControlPoint& base_cp,
                                              int maxDepth = 3) const;

  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      const Segment& base_segment, int maxDepth = 3) const;
  std::vector<image_t> GetNeighboringImageIds(const Segment& base_segment,
                                              int maxDepth = 3) const;

  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      const std::vector<Node>& base_nodes, int maxDepth = 3) const;
  std::vector<image_t> GetNeighboringImageIds(
      const std::vector<Node>& base_nodes, int maxDepth = 3) const;

  // utilities
  void AddControlPoint(const ControlPoint& cp);
  bool HasControlPoint(const ControlPoint& cp) const;
  void AddSegment(const Segment& segment);
  bool HasSegment(const Segment& segment) const;
  int NumControlPoints() const;
  int NumSegments() const;
  int NumNodes() const;

  void AddEdge(const ControlPoint& cp, const Segment& segment);
  void AddEdge(const ControlPoint& cp1, const ControlPoint& cp2);
  void AddEdge(const Segment& segment1, const Segment& segment2);
  int NumEdges() const;

  // tmp for pybind
  Node GetNode(const ControlPoint& cp) const;
  Node GetNode(const Segment& segment) const;

  std::map<int, ControlPointSequence> sequences;

 private:
  void AddNode(const Node& node);
  bool HasNode(const Node& node) const;
  void AddEdge(const Node& node1, const Node& node2);
  Node GetNode(const image_t image_id) const;
  void GetNeighboringNodes(const Node& node,
                           std::vector<Node>* neighbors) const;
  // customized bfs. max_depth limits the maximum control points that can be
  // traversed
  std::map<int, std::pair<timestamp_t, timestamp_t>> GetNeighboringRanges(
      std::queue<std::pair<Node, int>>& q,
      std::set<Node>& visited,
      int maxDepth = 3) const;

  // graph
  std::map<Node, std::set<Node>> g_nodes_;
  int num_edges_ = 0;

  // To help match cp with the same names
  std::map<std::string, std::vector<Node>> cp_name_to_nodes_;

  // image_id -> node
  std::map<image_t, Node> m_image_id_to_node_;
};

}  // namespace cp_partitioning

}  // namespace colmap
