#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <cstdint>
#include <sstream>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include "algorithms.h"
#include "graph.h"
#include "heap.h"

using FullEdges = std::vector<FullEdge>;
using ArgsType = std::vector<std::string>;
using EdgeID = std::uint_fast64_t;  // overwriting original definition

template <typename T>
class RandIntGen {
public:
    RandIntGen() = delete;
    RandIntGen(T n) : dis(0, n), gen(rd()) {}
    int operator()() {
        return dis(gen);
    }
private:
    std::random_device rd;
    std::uniform_int_distribution<T> dis;
    std::mt19937 gen;
};

// This function is to convert to template the function std::stoi
template <typename T>
T convertStringToInt(const std::string& str) {
    T result;
    std::stringstream ss(str);
    ss >> result;
    return result;
}

// Hash function
struct hashFunction {
  size_t operator()(const std::pair<NodeID,
                    NodeID> &x) const
  {
    return x.first ^ x.second;
  }
};

template <typename T>
void shuffle_vec(T it_begin, T it_end) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(it_begin, it_end, g);
}

class EdgeSet {
    public:
    EdgeSet(NodeID number_of_nodes, const FullEdges& edges, bool verbose=false) : n(number_of_nodes), nodes_sets(number_of_nodes) {
        auto progressBar = [](double progress, unsigned int barLength) {
            int pos = barLength * progress;
            std::clog << "\r[" << std::string(pos, '=') << std::string(barLength - pos, ' ') << "] " << static_cast<int>(progress * 100.0) << " %" << std::flush;
        };
        auto max_edges = edges.size();
        FullEdges::size_type counter = 0;
        for (auto const& [source, target, weight] : edges) {
            nodes_sets[source].insert(target);
            if (verbose) {
                progressBar(static_cast<double>(++counter) / max_edges, 50);
            }
        }
        if (verbose) std::clog << std::endl;
    }

    bool contains(const std::pair<NodeID, NodeID> &edge) const {
        if (edge.first >= n || edge.second >= n) return false;
        return nodes_sets[edge.first].find(edge.second) != nodes_sets[edge.first].end();
    }

    void insert(const std::pair<NodeID, NodeID> &edge) {
        if (edge.first >= n || edge.second >= n) return;
        nodes_sets[edge.first].insert(edge.second);
    }


    private:
    std::vector<std::unordered_set<NodeID>> nodes_sets;
    NodeID n;

};

FullEdges augment_graph(NodeID n, const FullEdges& edges, double p, Distance w) {
    EdgeID m_old_edges = edges.size();
    EdgeID max_new_edges = p <= 1.0 ? std::max(p * n * (n - 1) - m_old_edges, 0.0) : std::min(static_cast<EdgeID>(p * m_old_edges), static_cast<EdgeID>(n) * (n - 1)) - m_old_edges;

    std::clog << "Augmenting edges" << std::endl;

    auto progressBar = [](double progress, unsigned int barLength) {
        int pos = barLength * progress;
        std::clog << "\r[" << std::string(pos, '=') << std::string(barLength - pos, ' ') << "] " << static_cast<int>(progress * 100.0) << " %" << std::flush;
    };

    // std::clog << "Creating set of edges" << std::endl;
    // std::unordered_set<std::pair<NodeID, NodeID>, hashFunction> edges_set;
    // std::unordered_set<std::pair<NodeID, NodeID>, hashFunction>::size_type edge_count = 0;
    // edges_set.reserve(m_old_edges);
    // for (const auto& e : edges) {
    //     edges_set.insert({std::get<0>(e), std::get<1>(e)});
    //     progressBar(++edge_count / static_cast<double>(m_old_edges), 50);
    // }
    // std::clog << std::endl;

    EdgeSet edges_set(n, edges, true);

    // Deciding whether drawing the edges one by one or computing all the edges
    // first and then selecting a subset

    // Number or attempts in case we draw the edges one by one
    double exp_attempts = 0;
    for (EdgeID i = m_old_edges; i < n; i++) {
        exp_attempts += n / (n - i);
    }

    double all_possible_edges = static_cast<double>(n) * (n - 1);

    // Drawing one by one
    if (exp_attempts < all_possible_edges) {
        std::clog << "Drawing edges one by one" << std::endl;
        FullEdges new_edges;
        FullEdges::size_type new_edges_count = 0;
        RandIntGen<NodeID> gen(n - 1);
        while (new_edges.size() < max_new_edges) {
            NodeID i = gen();
            NodeID j = gen();
            if (i != j && !edges_set.contains({i, j})) {
            // if (i != j && edges_set.count({i, j}) == 0) {
                new_edges.emplace_back(i, j, w);
                edges_set.insert({i, j});
                progressBar(static_cast<double>(++new_edges_count) / max_new_edges, 50);
            }
        }
        std::clog << std::endl;
        return new_edges;
    }

    std::clog << "Computing all the edges" << std::endl;

    // Computing all the edges and then selecting a subset
    FullEdges new_edges;
    FullEdges::size_type new_edges_count = 0;
    for (NodeID i = 0; i < n; i++) {
        for (NodeID j = 0; j < n; j++) {
            if (i != j && !edges_set.contains({i, j})) {
            // if (i != j and edges_set.count({i, j}) == 0) {
                new_edges.push_back({i, j, w});
                progressBar(static_cast<double>(++new_edges_count) / max_new_edges, 50);
            }
        }
    }
    std::clog << std::endl;

    shuffle_vec(new_edges.begin(), new_edges.end());

    FullEdges all_edges(edges);
    all_edges.insert(all_edges.end(), new_edges.begin(), new_edges.begin() + static_cast<EdgeID>(max_new_edges));
    return all_edges;
}

std::vector<int> longest_paths(const Graph &graph) {
    NodeID n = graph.numberOfNodes();
    AddressableKHeap<4, NodeID, Distance> q(n);
    Distances d(n, c::infty);
    std::vector<int> count(n, 0);
    d[0] = 0;
    q.insert(0, 0);
    while (!q.empty()) {
        Distance dist;
        NodeID from;
        q.deleteMin(from, dist);
        if (dist > d[from]) continue;
        for (auto const &edge: graph.getEdgesOf(from)) {
            auto tentative_dist = d[from] + edge.weight;
            if (tentative_dist < d[edge.target]) {
                count[edge.target] = count[from] + 1;
                d[edge.target] = tentative_dist;
                q.insert(edge.target, tentative_dist);
            }
            else if (tentative_dist == d[edge.target])
                count[edge.target] = std::min(count[edge.target], count[from] + 1);
        }
    }

    return count;
}

void potential_transformation(Graph &graph) {
    std::cout << "graph.minWeight() = " << graph.minWeight() << std::endl;
    std::cout << "graph.maxWeight() = " << graph.maxWeight() << std::endl;
    long long int MAXX = std::max(graph.maxWeight(), -graph.minWeight());
    std::cout << "Applying potential transformation with MAX_POTENTIAL = " << MAXX << std::endl;

    long long int MIN_POTENTIAL = -MAXX, MAX_POTENTIAL = MAXX;
    RandIntGen<Distance> gen(MAX_POTENTIAL - MIN_POTENTIAL);
    NodeID n = graph.numberOfNodes();

    Distances potential(n, 0);
    for (NodeID i = 0; i < n; i++)
        potential[i] = MIN_POTENTIAL + gen();

    graph.applyPotential(potential);
}

std::pair<NodeID, FullEdges> get_edges(const Graph &graph) {
    NodeID n = graph.numberOfNodes();
    FullEdges edges;
    for (NodeID v = 0; v < n; v++) {
        for (auto const& e : graph.getEdgesOf(v)) {
            edges.emplace_back(v, e.target, e.weight);
        }
    }
    return {n, edges};
}

void graph_info(Graph &graph) {
    NodeID n = graph.numberOfNodes();
    EdgeID m = graph.numberOfEdges();
    
    // Count negative edges
    EdgeID negative_edges = 0;
    for (NodeID v = 0; v < n; v++) {
        for (auto const& edge : graph.getEdgesOf(v)) {
            if (edge.weight < 0) {
                negative_edges++;
            }
        }
    }
    
    // Print max edge weight and min edge weight
    std::cout << "Graph Info:" << std::endl;
    std::cout << "  Nodes: " << n << std::endl;
    std::cout << "  Edges: " << m << std::endl;
    std::cout << "  Negative Edges: " << negative_edges << std::endl;
    std::cout << "  Min Weight: " << graph.minWeight() << std::endl;
    std::cout << "  Max Weight: " << graph.maxWeight() << std::endl;
}

std::string appendToFilename(const std::string& input_name, const std::string& suffix) {
    std::string out_name = input_name;
    size_t slash_pos = out_name.find_last_of("/\\");
    size_t dot_pos = out_name.find_last_of('.');
    if (dot_pos != std::string::npos && (slash_pos == std::string::npos || dot_pos > slash_pos)) {
        out_name.insert(dot_pos, suffix);
    } else {
        out_name += suffix;
    }
    return out_name;
}

enum class OutputMode {
    MinusOne,
    Zero,
    TransformOnly
};

struct ModeConfig {
    OutputMode mode = OutputMode::MinusOne;
    std::string suffix = "_-1";
    bool add_edge = true;
    Distance edge_weight = -1;
};

ModeConfig parseMode(const std::string& mode_name) {
    if (mode_name == "-1") {
        return {OutputMode::MinusOne, "_-1", true, -1};
    }
    if (mode_name == "0") {
        return {OutputMode::Zero, "_0", true, 0};
    }
    if (mode_name == "t") {
        return {OutputMode::TransformOnly, "_t", false, 0};
    }

    throw std::invalid_argument("Invalid mode: expected -1, 0, or t");
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Invalid number of arguments: ./create_graph <graph type>" << std::endl;
        return 1;
    }

    ArgsType args;
    ModeConfig mode_config;


    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-m") {
            if (i + 1 >= argc) {
                std::cerr << "Invalid number of arguments for -m (need 1: -1, 0, or t)" << std::endl;
                return 1;
            }
            try {
                mode_config = parseMode(argv[i + 1]);
            } catch (const std::invalid_argument& ex) {
                std::cerr << ex.what() << std::endl;
                return 1;
            }
            i++;
            continue;
        }

        args.push_back(argv[i]);
    }

    for (auto it:args)
        std::cout << it << '\n';
    std::cout.flush();


    Graph graph =  readGraph(args[0]);

    auto opt_d = computeSSSP(SSSPAlg::LazyD, graph, 0);
    if (!opt_d.has_value())
        exit(-1);

    graph.applyPotential(opt_d.value());
    auto count = longest_paths(graph);

    NodeID n = graph.numberOfNodes();
    int maxx = -1;
    std::vector<NodeID> longest_path_nodes;
    for (NodeID i = 0; i < n; i++) {
        if (count[i] > maxx) {
            maxx = count[i];
            longest_path_nodes.clear();
            longest_path_nodes.push_back(i);
        }
        else if (count[i] == maxx) {
            longest_path_nodes.push_back(i);
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<std::size_t> dist(0, longest_path_nodes.size() - 1);
    NodeID ans = longest_path_nodes[dist(gen)];

    std::cout << "Longest path from node 0 is to node " << ans << " with length " << maxx << std::endl;
    std::cout<< "There are " << longest_path_nodes.size() << " nodes with the same longest path length" << std::endl;

    auto [n_edges, edges] = get_edges(graph);
    FullEdges output_edges = edges;
    if (mode_config.add_edge) {
        output_edges.emplace_back(ans, 0, mode_config.edge_weight);
    }

    Graph graph_w_cycle(n, output_edges);

    potential_transformation(graph_w_cycle);
    
    std::string out_name = appendToFilename(args[0], mode_config.suffix);
    std::ofstream fout(out_name);
    if (!fout) {
        std::cerr << "Could not open output file: " << out_name << std::endl;
        return 1;
    }
    graph_w_cycle.format_print(fout);

    graph_info(graph_w_cycle);
}