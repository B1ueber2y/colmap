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

#include "colmap/estimators/cp_partitioning.h"

#include "colmap/util/logging.h"

#include <queue>

namespace colmap {

namespace cp_partitioning {

ControlPointSequence::ControlPointSequence(
    const std::vector<ControlPoint>& control_points,
    const std::pair<timestamp_t, timestamp_t>& time_ranges)
    : control_points(control_points), time_ranges(time_ranges) {
  sequence_id = control_points[0].sequence_id;
  for (int i = 0; i < int(control_points.size()); ++i) {
    THROW_CHECK_EQ(control_points[i].cp_id, i);
    THROW_CHECK_EQ(control_points[i].sequence_id, sequence_id);
  }
  segments.push_back(Segment(sequence_id, 0, -1, 0));
  for (int i = 0; i < int(control_points.size()) - 1; ++i) {
    Segment segment(sequence_id, int(segments.size()), i, i + 1);
    segments.push_back(segment);
  }
  segments.push_back(Segment(sequence_id,
                             int(control_points.size()),
                             int(control_points.size()) - 1,
                             -1));
}

timestamp_t ControlPointSequence::GetSegmentStartTime(
    const int segment_id) const {
  Segment segment = segments[segment_id];
  if (segment.cp_id_left == -1)
    return time_ranges.first;
  else
    return control_points[segment.cp_id_left].timestamps.second;
}

timestamp_t ControlPointSequence::GetSegmentEndTime(
    const int segment_id) const {
  Segment segment = segments[segment_id];
  if (segment.cp_id_right == -1)
    return time_ranges.second;
  else
    return control_points[segment.cp_id_right].timestamps.first;
}

void ControlPointSequence::ImportImageTimestamps(
    const std::map<image_t, timestamp_t>& image_timestamps) {
  // set variables
  image_timestamps_ = image_timestamps;

  // generate indexing list
  std::map<timestamp_t, std::pair<NodeType, int>> start_times;
  for (auto& cp : control_points) {
    start_times.emplace(cp.timestamps.first,
                        std::make_pair(NodeType::CP, cp.cp_id));
  }
  for (auto& segment : segments) {
    timestamp_t starttime = GetSegmentStartTime(segment.segment_id);
    start_times.emplace(starttime,
                        std::make_pair(NodeType::SEGMENT, segment.segment_id));
  }

  // map images to cp / segment
  for (auto [image_id, timestamp] : image_timestamps) {
    THROW_CHECK_LE(time_ranges.first, timestamp);
    THROW_CHECK_GE(time_ranges.second, timestamp);
    auto it = start_times.upper_bound(timestamp);
    images_.emplace(image_id, it->second);
  }
}

std::pair<NodeType, int> ControlPointSequence::GetIndex(
    const image_t image_id) const {
  THROW_CHECK(images_.find(image_id) != images_.end());
  return images_.at(image_id);
}

std::vector<image_t> ControlPointSequence::GetImageIdsInsideTimeRange(
    const std::pair<timestamp_t, timestamp_t>& time_range) const {
  std::vector<image_t> image_ids;
  for (auto [image_id, timestamp] : image_timestamps_) {
    if (timestamp >= time_range.first && timestamp <= time_range.second)
      image_ids.push_back(image_id);
  }
  return image_ids;
}

void ControlPointSegmentGraph::AddNode(const Node& node) {
  if (!HasNode(node)) {
    g_nodes_.emplace(node, std::set<Node>());
  }
}

bool ControlPointSegmentGraph::HasNode(const Node& node) const {
  return g_nodes_.find(node) != g_nodes_.end();
}

Node ControlPointSegmentGraph::GetNode(const ControlPoint& cp) const {
  return Node(NodeType::CP, std::make_pair(cp.sequence_id, cp.cp_id));
}

Node ControlPointSegmentGraph::GetNode(const Segment& segment) const {
  return Node(NodeType::SEGMENT,
              std::make_pair(segment.sequence_id, segment.segment_id));
}

void ControlPointSegmentGraph::AddControlPoint(const ControlPoint& cp) {
  AddNode(GetNode(cp));
}

bool ControlPointSegmentGraph::HasControlPoint(const ControlPoint& cp) const {
  return HasNode(GetNode(cp));
}

void ControlPointSegmentGraph::AddSegment(const Segment& segment) {
  AddNode(GetNode(segment));
}

bool ControlPointSegmentGraph::HasSegment(const Segment& segment) const {
  return HasNode(GetNode(segment));
}

void ControlPointSegmentGraph::AddEdge(const Node& node1, const Node& node2) {
  THROW_CHECK(HasNode(node1));
  THROW_CHECK(HasNode(node2));
  g_nodes_.at(node1).insert(node2);
  g_nodes_.at(node2).insert(node1);
}

void ControlPointSegmentGraph::AddEdge(const ControlPoint& cp,
                                       const Segment& segment) {
  AddEdge(GetNode(cp), GetNode(segment));
}

void ControlPointSegmentGraph::AddEdge(const ControlPoint& cp1,
                                       const ControlPoint& cp2) {
  AddEdge(GetNode(cp1), GetNode(cp2));
}

void ControlPointSegmentGraph::AddEdge(const Segment& segment1,
                                       const Segment& segment2) {
  AddEdge(GetNode(segment1), GetNode(segment2));
}

void ControlPointSegmentGraph::GetNeighboringNodes(
    const Node& node, std::vector<Node>* neighbors) const {
  THROW_CHECK(HasNode(node));
  neighbors->clear();
  for (Node neighbor : g_nodes_.at(node)) {
    neighbors->push_back(neighbor);
  }
}

void ControlPointSegmentGraph::ImportSequence(
    const ControlPointSequence& sequence) {
  sequences.emplace(sequence.sequence_id, sequence);
  // add control points
  for (auto& cp : sequence.control_points) {
    AddControlPoint(cp);
    Node node_cp = GetNode(cp);
    if (cp_name_to_nodes_.find(cp.name) != cp_name_to_nodes_.end()) {
      // connect control points with the same name.
      for (Node& node_existing_cp : cp_name_to_nodes_.at(cp.name)) {
        AddEdge(node_cp, node_existing_cp);
      }
    } else {
      cp_name_to_nodes_.emplace(cp.name, std::vector<Node>());
    }
    cp_name_to_nodes_.at(cp.name).push_back(node_cp);
  }
  // add segments and connect with cp
  for (auto& segment : sequence.segments) {
    AddSegment(segment);
    if (segment.cp_id_left != -1) {
      AddEdge(sequence.control_points[segment.cp_id_left], segment);
    }
    if (segment.cp_id_right != -1) {
      AddEdge(sequence.control_points[segment.cp_id_right], segment);
    }
  }
}

void ControlPointSegmentGraph::ImportSequenceMatching(
    const SequenceMatching& matches) {
  THROW_CHECK(sequences.find(matches.sequence_id_1) != sequences.end());
  ControlPointSequence& seq1 = sequences.at(matches.sequence_id_1);
  THROW_CHECK(sequences.find(matches.sequence_id_2) != sequences.end());
  ControlPointSequence& seq2 = sequences.at(matches.sequence_id_2);

  for (auto& match : matches.matches) {
    auto index1 = seq1.GetIndex(match.first);
    Node node1 =
        Node(index1.first, std::make_pair(seq1.sequence_id, index1.second));
    auto index2 = seq2.GetIndex(match.second);
    Node node2 =
        Node(index2.first, std::make_pair(seq2.sequence_id, index2.second));
    AddEdge(node1, node2);
  }
}

std::map<int, std::pair<timestamp_t, timestamp_t>>
ControlPointSegmentGraph::GetNeighboringRanges(
    std::queue<std::pair<Node, int>>& q,
    std::set<Node>& visited,
    int maxDepth) const {
  // customized bfs
  while (!q.empty()) {
    auto [node, depth] = q.front();
    q.pop();
    if (depth >= maxDepth) {
      continue;
    }
    // traverse
    std::vector<Node> neighbors;
    GetNeighboringNodes(node, &neighbors);
    for (Node neighbor : neighbors) {
      if (visited.find(neighbor) != visited.end()) continue;
      visited.insert(neighbor);
      // increment depth only for segment - cp edge
      if (node.first == NodeType::SEGMENT && neighbor.first == NodeType::CP) {
        q.push({neighbor, depth + 1});
      }
      // otherwise keep the depth unchanged
      else {
        q.push({neighbor, depth});
      }
    }
  }

  // Get result
  std::map<int, std::pair<timestamp_t, timestamp_t>> res;
  for (auto& node : visited) {
    // get timestamps
    timestamp_t timestamp_low, timestamp_high;
    int sequence_id = node.second.first;
    if (node.first == NodeType::CP) {
      int cp_id = node.second.second;
      timestamp_low =
          sequences.at(sequence_id).control_points[cp_id].timestamps.first;
      timestamp_high =
          sequences.at(sequence_id).control_points[cp_id].timestamps.second;
    } else {
      int segment_id = node.second.second;
      timestamp_low = sequences.at(sequence_id).GetSegmentStartTime(segment_id);
      timestamp_high = sequences.at(sequence_id).GetSegmentEndTime(segment_id);
    }
    if (res.find(sequence_id) == res.end()) {
      res.emplace(sequence_id, std::make_pair(timestamp_low, timestamp_high));
    } else {
      if (res[sequence_id].first > timestamp_low) {
        res[sequence_id].first = timestamp_low;
      }
      if (res[sequence_id].first < timestamp_high) {
        res[sequence_id].second = timestamp_high;
      }
    }
  }
  return res;
}

std::map<int, std::pair<timestamp_t, timestamp_t>>
ControlPointSegmentGraph::GetNeighboringRanges(const ControlPoint& base_cp,
                                               int maxDepth) const {
  Node node = GetNode(base_cp);
  std::queue<std::pair<Node, int>> q;
  std::set<Node> visited;
  visited.insert(node);
  q.push({node, 0});
  return GetNeighboringRanges(q, visited, maxDepth);
}

std::map<int, std::pair<timestamp_t, timestamp_t>>
ControlPointSegmentGraph::GetNeighboringRanges(const Segment& base_segment,
                                               int maxDepth) const {
  std::queue<std::pair<Node, int>> q;
  std::set<Node> visited;
  Node node = GetNode(base_segment);
  visited.insert(node);
  q.push({node, 0});
  return GetNeighboringRanges(q, visited, maxDepth);
}

std::vector<image_t> ControlPointSegmentGraph::GetNeighboringImageIds(
    const ControlPoint& base_cp, int maxDepth) const {
  std::vector<image_t> image_ids;
  std::map<int, std::pair<timestamp_t, timestamp_t>> time_ranges =
      GetNeighboringRanges(base_cp, maxDepth);
  for (auto& [seq_id, time_range] : time_ranges) {
    std::vector<image_t> image_ids_in_seq =
        sequences.at(seq_id).GetImageIdsInsideTimeRange(time_range);
    image_ids.insert(
        image_ids.end(), image_ids_in_seq.begin(), image_ids_in_seq.end());
  }
  return image_ids;
}

std::vector<image_t> ControlPointSegmentGraph::GetNeighboringImageIds(
    const Segment& base_segment, int maxDepth) const {
  std::vector<image_t> image_ids;
  std::map<int, std::pair<timestamp_t, timestamp_t>> time_ranges =
      GetNeighboringRanges(base_segment, maxDepth);
  for (auto& [seq_id, time_range] : time_ranges) {
    std::vector<image_t> image_ids_in_seq =
        sequences.at(seq_id).GetImageIdsInsideTimeRange(time_range);
    image_ids.insert(
        image_ids.end(), image_ids_in_seq.begin(), image_ids_in_seq.end());
  }
  return image_ids;
}

}  // namespace cp_partitioning

}  // namespace colmap