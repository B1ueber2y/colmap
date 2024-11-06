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

ControlPointIndex ControlPointSegmentGraph::GetControlPointIndex(
    const ControlPoint& cp) const {
  return ControlPointIndex(cp.sequence_id, cp.cp_id);
}

ControlPoint ControlPointSegmentGraph::GetControlPoint(
    const ControlPointIndex& cp_index) const {
  int sequence_id = cp_index.first;
  int cp_id = cp_index.second;
  THROW_CHECK(sequences_.find(sequence_id) != sequences_.end());
  return sequences_.at(sequence_id)->control_points[cp_id];
}

SegmentIndex ControlPointSegmentGraph::GetSegmentIndex(
    const Segment& segment) const {
  return SegmentIndex(segment.sequence_id, segment.segment_id);
}

Segment ControlPointSegmentGraph::GetSegment(
    const SegmentIndex& segment_index) const {
  int sequence_id = segment_index.first;
  int segment_id = segment_index.second;
  THROW_CHECK(sequences_.find(sequence_id) != sequences_.end());
  return sequences_.at(sequence_id)->segments[segment_id];
}

void ControlPointSegmentGraph::AddControlPoint(const ControlPoint& cp) {
  auto cp_index = GetControlPointIndex(cp);
  m_cp_to_cp.emplace(cp_index, std::set<ControlPointIndex>());
  m_cp_to_segment.emplace(cp_index, std::set<ControlPointIndex>());
}

bool ControlPointSegmentGraph::HasControlPoint(const ControlPoint& cp) const {
  auto cp_index = GetControlPointIndex(cp);
  return HasControlPoint(cp_index);
}

bool ControlPointSegmentGraph::HasControlPoint(
    const ControlPointIndex& cp_index) const {
  return m_cp_to_cp.find(cp_index) != m_cp_to_cp.end();
}

void ControlPointSegmentGraph::AddSegment(const Segment& segment) {
  auto segment_index = GetSegmentIndex(segment);
  m_segment_to_cp.emplace(segment_index, std::set<ControlPointIndex>());
  m_segment_to_segment.emplace(segment_index, std::set<SegmentIndex>());
}

bool ControlPointSegmentGraph::HasSegment(const Segment& segment) const {
  auto segment_index = GetSegmentIndex(segment);
  return HasSegment(segment_index);
}

bool ControlPointSegmentGraph::HasSegment(
    const SegmentIndex& segment_index) const {
  return m_segment_to_segment.find(segment_index) != m_segment_to_segment.end();
}

void ControlPointSegmentGraph::AddEdge(const ControlPoint& cp,
                                       const Segment& segment) {
  THROW_CHECK(HasControlPoint(cp));
  auto cp_index = GetControlPointIndex(cp);
  THROW_CHECK(HasSegment(segment));
  auto segment_index = GetSegmentIndex(segment);
  m_cp_to_segment.at(cp_index).insert(segment_index);
  m_segment_to_cp.at(segment_index).insert(cp_index);
}

void ControlPointSegmentGraph::AddEdge(const ControlPoint& cp1,
                                       const ControlPoint& cp2) {
  THROW_CHECK(HasControlPoint(cp1));
  auto cp1_index = GetControlPointIndex(cp1);
  THROW_CHECK(HasControlPoint(cp2));
  auto cp2_index = GetControlPointIndex(cp2);
  THROW_CHECK(cp1_index != cp2_index);
  m_cp_to_cp.at(cp1_index).insert(cp2_index);
  m_cp_to_cp.at(cp2_index).insert(cp1_index);
}

void ControlPointSegmentGraph::AddEdge(const Segment& segment1,
                                       const Segment& segment2) {
  THROW_CHECK(HasSegment(segment1));
  auto segment1_index = GetSegmentIndex(segment1);
  THROW_CHECK(HasSegment(segment2));
  auto segment2_index = GetSegmentIndex(segment2);
  THROW_CHECK(segment1_index != segment2_index);
  m_segment_to_segment.at(segment1_index).insert(segment2_index);
  m_segment_to_segment.at(segment2_index).insert(segment1_index);
}

void ControlPointSegmentGraph::GetNeighboringSegmentsFromCP(
    const ControlPointIndex& cp_index,
    std::vector<SegmentIndex>* neighboring_segments) const {
  THROW_CHECK(HasControlPoint(cp_index));
  neighboring_segments->clear();
  for (auto neighbor : m_cp_to_segment.at(cp_index)) {
    neighboring_segments->push_back(neighbor);
  }
}

void ControlPointSegmentGraph::GetNeighboringSegmentsFromSegment(
    const SegmentIndex& segment_index,
    std::vector<SegmentIndex>* neighboring_segments) const {
  THROW_CHECK(HasSegment(segment_index));
  neighboring_segments->clear();
  for (auto neighbor : m_segment_to_segment.at(segment_index)) {
    neighboring_segments->push_back(neighbor);
  }
}

void ControlPointSegmentGraph::GetNeighboringControlPointsFromCP(
    const ControlPointIndex& cp_index,
    std::vector<ControlPointIndex>* neighboring_cps) const {
  THROW_CHECK(HasControlPoint(cp_index));
  neighboring_cps->clear();
  for (auto neighbor : m_cp_to_cp.at(cp_index)) {
    neighboring_cps->push_back(neighbor);
  }
}

void ControlPointSegmentGraph::GetNeighboringControlPointsFromSegment(
    const SegmentIndex& segment_index,
    std::vector<ControlPointIndex>* neighboring_cps) const {
  THROW_CHECK(HasSegment(segment_index));
  neighboring_cps->clear();
  for (auto neighbor : m_segment_to_cp.at(segment_index)) {
    neighboring_cps->push_back(neighbor);
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
ControlPointSegmentGraph::GetNeighboringRanges(const ControlPoint& base_cp,
                                               int maxDepth) const {
  // BFS
  auto base_cp_index = GetControlPointIndex(base_cp);
  std::queue<std::pair<ControlPointIndex, int>> q;
  std::set<ControlPointIndex> visited;
  visited.insert(base_cp_index);
  q.push({base_cp_index, 0});

  while (!q.empty()) {
    auto [node, depth] = q.front();
    q.pop();
    if (depth >= maxDepth) {
      continue;
    }
    // traverse
    std::vector<SegmentIndex> segments;
    GetNeighboringSegmentsFromCP(node, &segments);
    for (SegmentIndex segment_index : segments) {
      std::vector<ControlPointIndex> neighbors;
      GetNeighboringControlPointsFromSegment(segment_index, &neighbors);
      for (ControlPointIndex neighbor : neighbors) {
        if (visited.find(neighbor) != visited.end()) continue;
        visited.insert(neighbor);
        q.push({neighbor, depth + 1});
      }
    }
  }

  // Get result
  std::map<int, std::pair<timestamp_t, timestamp_t>> res;
  for (auto& node : visited) {
    auto cp = GetControlPoint(node);
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
ControlPointSegmentGraph::GetNeighboringRanges(const Segment& base_segment,
                                               int maxDepth) const {
  // BFS
  std::queue<std::pair<ControlPointIndex, int>> q;
  std::set<ControlPointIndex> visited;
  std::vector<ControlPointIndex> cps;
  GetNeighboringControlPointsFromSegment(GetSegmentIndex(base_segment), &cps);
  for (ControlPointIndex cp_index : cps) {
    visited.insert(cp_index);
    q.push({cp_index, 0});
  }

  while (!q.empty()) {
    auto [node, depth] = q.front();
    q.pop();
    if (depth >= maxDepth) {
      continue;
    }
    // traverse
    std::vector<SegmentIndex> segments;
    GetNeighboringSegmentsFromCP(node, &segments);
    for (SegmentIndex segment_index : segments) {
      std::vector<ControlPointIndex> neighbors;
      GetNeighboringControlPointsFromSegment(segment_index, &neighbors);
      for (ControlPointIndex neighbor : neighbors) {
        if (visited.find(neighbor) != visited.end()) continue;
        visited.insert(neighbor);
        q.push({neighbor, depth + 1});
      }
    }
  }

  // Get result
  std::map<int, std::pair<timestamp_t, timestamp_t>> res;
  for (auto& node : visited) {
    auto cp = GetControlPoint(node);
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

}  // namespace colmap
