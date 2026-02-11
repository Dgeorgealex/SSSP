//
// Created by adumi on 1/27/26.
//

#ifndef NEGATIVEWEIGHTSHORTESTPATH_PAD_H
#define NEGATIVEWEIGHTSHORTESTPATH_PAD_H

#include <optional>
#include "defs.h"
#include "graph.h"
#include "heap.h"

namespace pad {
    using GraphHeap = AddressableKHeap<4, NodeID, Distance>;

    class PADAlg {
    public:
        PADAlg() = default;

        std::optional<Distances> runMainAlg(Graph &graph, Distance diameter, int level = 0);
    };
} // pad

EdgeID grow_ball(Distance new_d, int &start, int end, const std::vector<Distance> &d,
                 const std::vector<NodeID> &order,
                 const std::vector<EdgeID> &vol, Distance diameter);

std::vector<Graph> padded_decomposition(Graph &graph, Distance diameter);

std::pair<int, int> grow_ball_heavy(const Graph &graph, Distance diameter, std::vector<NodeID> &order, NodeID s,
                                    Orientation orientation);

std::vector<Graph> padded_decomposition_heavy(const Graph &graph, Distance diameter, NodeID s);
#endif //NEGATIVEWEIGHTSHORTESTPATH_PAD_H
