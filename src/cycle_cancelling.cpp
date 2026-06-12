#include "algorithms.h"
#include "gor.h"
#include "pad.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <ranges>
#include <string>
#include <vector>

#include "measurement_tool.h"

void read_graph(std::string filename);

void init_feasible();

Graph get_residual_graph();

std::optional<std::vector<NodeID> > get_cycle(Graph &graph, SSSPAlg alg);

void cancel_cycle(std::vector<NodeID> cycle);

struct OrigEdge {
    NodeID source;
    NodeID target;
    Distance cost;
    Distance flow;
};

struct FlowNetwork {
    NodeID number_of_nodes = 0;
    std::vector<Distance> balance;
    std::vector<OrigEdge> edges;
    std::vector<std::vector<EdgeID> > out_edges;
};

FlowNetwork g_network;

EdgeID find_edge_id(NodeID u, NodeID v) {
    for (EdgeID id: g_network.out_edges[u]) {
        if (g_network.edges[id].target == v) {
            return id;
        }
    }
    return static_cast<EdgeID>(-1);
}

int iterations;
Distance old_cost = c::infty;

void quick_correctness_check() {
    std::vector<Distance> check(g_network.number_of_nodes, 0);
    Distance cost = 0;
    for (auto it: g_network.edges) {
        check[it.source] -= it.flow;
        check[it.target] += it.flow;
        cost += it.flow;
    }

    for (int i = 0; i < g_network.number_of_nodes; ++i)
        if (check[i] != g_network.balance[i]) {
            std::cout << "Bad relaxation\n";
            std::cout.flush();
            exit(-1);
        }

    if (cost >= old_cost) {
        std::cout << "Bad cost\n";
        std::cout.flush();
        exit(-1);
    }

    old_cost = cost;
}

int main() {
    // First incremental step: parse and validate the instance format.
    read_graph("../data/Pictures/graph.txt");
    std::cout << "Loaded graph with " << g_network.number_of_nodes << " nodes and "
            << g_network.edges.size() << " edges.\n";

    init_feasible();

    quick_correctness_check();

    while (true) {
        iterations++;
        auto graph = get_residual_graph();
        auto cycle = get_cycle(graph, SSSPAlg::PADSCALING);

        if (!cycle.has_value())
            break;

        cancel_cycle(cycle.value());
    }
    return 0;
}


void read_graph(std::string filename) {
    std::ifstream file(filename);
    NodeID n;
    EdgeID m;
    file >> n >> m;

    g_network.number_of_nodes = n;
    g_network.balance.resize(n, 0);
    g_network.out_edges.assign(n, {});
    g_network.edges.reserve(m);

    for (NodeID v = 0; v < n; ++v) {
        file >> g_network.balance[v];
    }

    for (EdgeID i = 0; i < m; ++i) {
        NodeID u;
        NodeID v;
        file >> u >> v;

        const EdgeID edge_id = g_network.edges.size();
        g_network.edges.push_back(OrigEdge{u, v, 1, 0});
        g_network.out_edges[u].push_back(edge_id);
    }
}

void init_feasible() {
    const NodeID n = g_network.number_of_nodes;

    for (auto &edge: g_network.edges) {
        edge.flow = 0;
    }

    if (n == 0) return;

    const NodeID root = 0;

    // Undirected adjacency (assume for every u->v there exists v->u)
    std::vector<std::vector<NodeID> > adj(n);
    for (const auto &e: g_network.edges) {
        adj[e.source].push_back(e.target);
        adj[e.target].push_back(e.source);
    }

    std::vector<NodeID> parent(n, n);
    std::vector<std::vector<NodeID> > children(n);
    std::vector<NodeID> order;
    std::vector<bool> visited(n, false);

    std::queue<NodeID> q;
    q.push(root);
    visited[root] = true;
    while (!q.empty()) {
        NodeID v = q.front();
        q.pop();
        order.push_back(v);
        for (NodeID to: adj[v]) {
            if (visited[to]) continue;
            visited[to] = true;
            parent[to] = v;
            children[v].push_back(to);
            q.push(to);
        }
    }

    std::vector<Distance> current_balance = g_network.balance;
    for (unsigned long v: std::views::reverse(order)) {
        for (auto c: children[v]) {
            if (current_balance[c] < 0) {
                EdgeID id = find_edge_id(c, v);
                g_network.edges[id].flow = -current_balance[c];
            } else if (current_balance[c] > 0) {
                EdgeID id = find_edge_id(v, c);
                g_network.edges[id].flow = current_balance[c];
            }
            current_balance[v] += current_balance[c];
            current_balance[c] = 0;
        }
    }
}

Graph get_residual_graph() {
    // Build a residual graph: forward arcs (u->v) with cost = +orig.cost
    // and reverse arcs (v->u) with cost = -orig.cost if flow > 0.
    const NodeID n = g_network.number_of_nodes;
    using FullEdge = std::tuple<NodeID, NodeID, Distance>;
    std::vector<FullEdge> edges;
    edges.reserve(g_network.edges.size() * 2);

    for (EdgeID id = 0; id < g_network.edges.size(); ++id) {
        const OrigEdge &e = g_network.edges[id];
        // forward residual
        edges.emplace_back(e.source, e.target, e.cost);
        // reverse residual exists only if we have positive flow
        if (e.flow > 0) {
            edges.emplace_back(e.target, e.source, -e.cost);
        }
    }

    return Graph(n, edges);
}

std::optional<std::vector<NodeID> > get_cycle(Graph &graph, SSSPAlg alg) {
    MEASUREMENT::start(EXP::INNER_LOOP_ALL);
    auto result1 = PADSCALING(graph, 0);
    MEASUREMENT::stop(EXP::INNER_LOOP_ALL);
    std::cout << "PAD: \n";
    MEASUREMENT::print(EXP::INNER_LOOP_ALL);
    MEASUREMENT::reset(EXP::INNER_LOOP_ALL);

    // MEASUREMENT::start(EXP::INNER_LOOP_ALL);
    // auto result2 = gor_with_cycles(graph, static_cast<NodeID>(0));
    // MEASUREMENT::stop(EXP::INNER_LOOP_ALL);
    // std::cout << "GOR: \n";
    // MEASUREMENT::print(EXP::INNER_LOOP_ALL);
    // MEASUREMENT::reset(EXP::INNER_LOOP_ALL);

    if (std::holds_alternative<Distances>(result1))
        return {};

    // if (std::holds_alternative<Distances>(result1))
    //     exit(-1);

    return std::get<std::vector<NodeID> >(result1);
}

void cancel_cycle(std::vector<NodeID> cycle) {
    if (cycle.empty()) return;

    const size_t k = cycle.size();
    const Distance INF = std::numeric_limits<Distance>::max() / 4;

    // For each residual arc in the cycle store (orig_edge_id, is_forward)
    std::vector<std::pair<EdgeID, bool> > arc_info(k);
    Distance bottleneck = INF;

    for (size_t i = 0; i < k; ++i) {
        NodeID u = cycle[i];
        NodeID v = cycle[(i + 1) % k];

        EdgeID fwd = find_edge_id(u, v);
        EdgeID rev = find_edge_id(v, u);

        if (rev != static_cast<EdgeID>(-1) && g_network.edges[rev].flow > 0) {
            arc_info[i] = {rev, false};
            bottleneck = std::min(bottleneck, g_network.edges[rev].flow);
        } else if (fwd != static_cast<EdgeID>(-1)) {
            arc_info[i] = {fwd, true};
        } else {
            // malformed cycle w.r.t. original graph
            std::cout << "Malformed cycle";
            std::cout.flush();
            exit(-1);
        }
    }

    // Apply augmentation of size `bottleneck`
    int better = 0;
    for (size_t i = 0; i < k; ++i) {
        EdgeID id = arc_info[i].first;
        bool is_fwd = arc_info[i].second;
        if (is_fwd) {
            g_network.edges[id].flow += bottleneck;
            better -= 1;
        } else {
            // reverse -> reduce flow on original edge
            if (g_network.edges[id].flow <= bottleneck) g_network.edges[id].flow = 0;
            else g_network.edges[id].flow -= bottleneck;

            better += 1;
        }
    }

    quick_correctness_check();

    std::cout << iterations << ' ' << cycle.size() << ' ' << better << ' ' << bottleneck << '\n';
}
