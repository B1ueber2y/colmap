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
    const std::vector<ControlPoint>& control_points)
    : control_points(control_points) {
  sequence_id = control_points[0].sequence_id;
  for (int i = 0; i < int(control_points.size()); ++i) {
    THROW_CHECK_EQ(control_points[i].cp_id, i);
    THROW_CHECK_EQ(control_points[i].sequence_id, sequence_id);
  }
  for (int i = 0; i < int(control_points.size()) - 1; ++i) {
    Segment segment(sequence_id, i, i, i + 1);
    segments.push_back(segment);
  }
}

void ControlPointSequence::ImportImages(
    const std::map<image_t, timestamp_t>& images) {
  // TODO: impl
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

void ControlPointSegmentGraph::ImportSequence(ControlPointSequence* sequence) {
  sequences_.emplace(sequence->sequence_id, sequence);
  for (auto& cp : sequence->control_points) {
    AddControlPoint(cp);
  }
  for (auto& segment : sequence->segments) {
    AddSegment(segment);
  }
  for (size_t i = 0; i < sequence->control_points.size() - 1; ++i) {
    AddEdge(sequence->control_points[i], sequence->segments[i]);
    AddEdge(sequence->control_points[i + 1], sequence->segments[i]);
  }
}

void ControlPointSegmentGraph::ImportSequenceMatching(
    const SequenceMatching& matches) {
  // TODO: impl later
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
    if (node.first != NodeType::CP) continue;
    int sequence_id = node.second.first;
    int cp_id = node.second.second;
    ControlPoint cp = sequences_.at(sequence_id)->control_points[cp_id];
    if (res.find(cp.sequence_id) == res.end()) {
      res.emplace(cp.sequence_id, cp.timestamps);
    } else {
      if (res[cp.sequence_id].first > cp.timestamps.first) {
        res[cp.sequence_id].first = cp.timestamps.first;
      }
      if (res[cp.sequence_id].first < cp.timestamps.first) {
        res[cp.sequence_id].second = cp.timestamps.second;
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

}  // namespace cp_partitioning

}  // namespace colmap
