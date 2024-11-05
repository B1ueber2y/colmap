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
    Segment segment(sequence_id, i, i + 1);
    segments.push_back(segment);
  }
}

void ControlPointSegmentGraph::AddControlPoint(ControlPoint cp) {
  m_cp_to_cp.emplace(cp, std::vector<ControlPoint*>());
  m_cp_to_segment.emplace(cp, std::vector<Segment*>());
}

bool ControlPointSegmentGraph::HasControlPoint(ControlPoint cp) const {
  return m_cp_to_cp.find(cp) != m_cp_to_cp.end();
}

void ControlPointSegmentGraph::AddSegment(Segment segment) {
  m_segment_to_cp.emplace(segment, std::vector<ControlPoint*>());
  m_segment_to_segment.emplace(segment, std::vector<Segment*>());
}

bool ControlPointSegmentGraph::HasSegment(Segment segment) const {
  return m_segment_to_segment.find(segment) != m_segment_to_segment.end();
}

void ControlPointSegmentGraph::AddEdge(ControlPoint* cp, Segment* segment) {
  THROW_CHECK(HasControlPoint(*cp));
  THROW_CHECK(HasSegment(*segment));
  m_cp_to_segment.at(*cp).push_back(segment);
  m_segment_to_cp.at(*segment).push_back(cp);
  // TODO: handle repetitive edges
}

void ControlPointSegmentGraph::AddEdge(ControlPoint* cp1, ControlPoint* cp2) {
  THROW_CHECK(HasControlPoint(*cp1));
  THROW_CHECK(HasControlPoint(*cp2));
  THROW_CHECK(cp1 != cp2);
  m_cp_to_cp.at(*cp1).push_back(cp2);
  m_cp_to_cp.at(*cp2).push_back(cp1);
}

void ControlPointSegmentGraph::AddEdge(Segment* segment1, Segment* segment2) {
  THROW_CHECK(HasSegment(*segment1));
  THROW_CHECK(HasSegment(*segment2));
  THROW_CHECK(segment1 != segment2);
  m_segment_to_segment.at(*segment1).push_back(segment2);
  m_segment_to_segment.at(*segment2).push_back(segment1);
}

void ControlPointSegmentGraph::GetNeighboringSegments(
    ControlPoint cp, std::vector<Segment*>* neighboring_segments) const {
  THROW_CHECK(HasControlPoint(cp));
  neighboring_segments->clear();
  for (Segment* seg : m_cp_to_segment.at(cp)) {
    neighboring_segments->push_back(seg);
  }
}

void ControlPointSegmentGraph::GetNeighboringSegments(
    Segment segment, std::vector<Segment*>* neighboring_segments) const {
  THROW_CHECK(HasSegment(segment));
  neighboring_segments->clear();
  for (Segment* seg : m_segment_to_segment.at(segment)) {
    neighboring_segments->push_back(seg);
  }
}

void ControlPointSegmentGraph::GetNeighboringControlPoints(
    ControlPoint cp, std::vector<ControlPoint*>* neighboring_cps) const {
  THROW_CHECK(HasControlPoint(cp));
  neighboring_cps->clear();
  for (ControlPoint* node : m_cp_to_cp.at(cp)) {
    neighboring_cps->push_back(node);
  }
}

void ControlPointSegmentGraph::GetNeighboringControlPoints(
    Segment segment, std::vector<ControlPoint*>* neighboring_cps) const {
  THROW_CHECK(HasSegment(segment));
  neighboring_cps->clear();
  for (ControlPoint* node : m_segment_to_cp.at(segment)) {
    neighboring_cps->push_back(node);
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
    AddEdge(&sequence->control_points[i], &sequence->segments[i]);
    AddEdge(&sequence->control_points[i + 1], &sequence->segments[i]);
  }
}

void ControlPointSegmentGraph::ImportCrossSequenceMatching(
    const CrossSequenceMatching& matches) {
  // TODO: impl later
}

std::map<int, std::pair<time_t, time_t>>
ControlPointSegmentGraph::GetNeighboringRanges(ControlPoint base_cp,
                                               int maxDepth) const {
  // BFS
  std::queue<std::pair<ControlPoint, int>> q;
  std::set<ControlPoint> visited;
  visited.insert(base_cp);
  q.push({base_cp, 0});

  while (!q.empty()) {
    auto [node, depth] = q.front();
    q.pop();
    if (depth >= maxDepth) {
      continue;
    }
    // traverse
    std::vector<Segment*> segments;
    GetNeighboringSegments(node, &segments);
    for (Segment* segment : segments) {
      std::vector<ControlPoint*> neighbors;
      GetNeighboringControlPoints(*segment, &neighbors);
      for (ControlPoint* neighbor : neighbors) {
        if (visited.find(*neighbor) != visited.end()) continue;
        visited.insert(*neighbor);
        q.push({*neighbor, depth + 1});
      }
    }
  }

  // Get result
  std::map<int, std::pair<time_t, time_t>> res;
  for (auto& node : visited) {
    if (res.find(node.sequence_id) == res.end()) {
      res.emplace(node.sequence_id, node.timestamps);
    } else {
      if (res[node.sequence_id].first > node.timestamps.first) {
        res[node.sequence_id].first = node.timestamps.first;
      }
      if (res[node.sequence_id].first < node.timestamps.first) {
        res[node.sequence_id].second = node.timestamps.second;
      }
    }
  }
  return res;
}

std::map<int, std::pair<time_t, time_t>>
ControlPointSegmentGraph::GetNeighboringRanges(Segment base_segment,
                                               int maxDepth) const {
  // BFS
  std::queue<std::pair<ControlPoint, int>> q;
  std::set<ControlPoint> visited;
  std::vector<ControlPoint*> cps;
  GetNeighboringControlPoints(base_segment, &cps);
  for (ControlPoint* cp : cps) {
    visited.insert(*cp);
    q.push({*cp, 0});
  }

  while (!q.empty()) {
    auto [node, depth] = q.front();
    q.pop();
    if (depth >= maxDepth) {
      continue;
    }
    // traverse
    std::vector<Segment*> segments;
    GetNeighboringSegments(node, &segments);
    for (Segment* segment : segments) {
      std::vector<ControlPoint*> neighbors;
      GetNeighboringControlPoints(*segment, &neighbors);
      for (ControlPoint* neighbor : neighbors) {
        if (visited.find(*neighbor) != visited.end()) continue;
        visited.insert(*neighbor);
        q.push({*neighbor, depth + 1});
      }
    }
  }

  // Get result
  std::map<int, std::pair<time_t, time_t>> res;
  for (auto& node : visited) {
    if (res.find(node.sequence_id) == res.end()) {
      res.emplace(node.sequence_id, node.timestamps);
    } else {
      if (res[node.sequence_id].first > node.timestamps.first) {
        res[node.sequence_id].first = node.timestamps.first;
      }
      if (res[node.sequence_id].first < node.timestamps.first) {
        res[node.sequence_id].second = node.timestamps.second;
      }
    }
  }
  return res;
}

}  // namespace colmap
