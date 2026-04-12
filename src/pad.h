//
// Created by adumi on 1/27/26.
//

#ifndef NEGATIVEWEIGHTSHORTESTPATH_PAD_H
#define NEGATIVEWEIGHTSHORTESTPATH_PAD_H

#include <optional>
#include "graph.h"
#include "heap.h"
#include <ranges>
#include "defs.h"
#include "bcf.h"
#include "algorithms.h"

namespace pad {
    using GraphHeap = AddressableKHeap<4, NodeID, Distance>;

    class PADAlg {
    public:
        PADAlg() = default;

        std::optional<Distances> runMainAlg(Graph &graph, Distance diameter, int level = 0);
    };

    bool fast_admissible_graph_check(const Graph &graph, const Distances &potential);

    std::optional<Distances> runLazyDijkstra(const Graph& graph, const Distances& potential, Distance diameter, int max_rounds);

    std::optional<Distances> scaling_early_finish(const Graph &graph, const Graph &current_graph, NodeID source);

    struct PADStats {
        bool in_padding = false;
        int scaling_iterations = 0;
        int max_recursion_level = 0;
        Distance final_minW = static_cast<Distance>(-1);
        int decomposition_calls = 0;
        int decomposition_calls_with_padding = 0;
        int heavy_calls = 0;

        void reset() {
            in_padding = false;
            scaling_iterations = 0;
            max_recursion_level = 0;
            final_minW = static_cast<Distance>(-1);
            decomposition_calls = 0;
            decomposition_calls_with_padding = 0;
            heavy_calls = 0;
        }
    };

    inline PADStats stats{};
} // pad

EdgeID grow_ball(Distance new_d, int &start, int end, const std::vector<Distance> &d,
                 const std::vector<NodeID> &order,
                 const std::vector<EdgeID> &vol, Distance diameter);

std::vector<Graph> padded_decomposition(Graph &graph, Distance diameter);

std::pair<int, int> grow_ball_heavy(const Graph &graph, Distance diameter, std::vector<NodeID> &order, NodeID s,
                                    Orientation orientation);

std::vector<Graph> padded_decomposition_heavy(const Graph &graph, Distance diameter, NodeID s);

void padding_check(const Graph &graph, const std::vector<bool> &u, const std::vector<bool> &pad, Distance diameter, Orientation orientation);
#endif //NEGATIVEWEIGHTSHORTESTPATH_PAD_H
