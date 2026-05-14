#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <cstdint>
#include <sstream>
#include <queue>
#include <numeric>
#include <cmath>
#include <fstream>

#include "graph.h"
#include "heap.h"

using FullEdges = std::vector<FullEdge>;
using GraphInfo = std::tuple<NodeID, FullEdges>;
using ArgsType = std::vector<std::string>;
using GraphCreatorType = GraphInfo (*)(const ArgsType&); // Define the function type
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

FullEdges getRandomSubset(const FullEdges& vec, NodeID k) {
    // Create a vector of indices from 0 to vec.size() - 1
    std::vector<NodeID> indices(vec.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Shuffle the indices randomly
    shuffle_vec(indices.begin(), indices.end());

    // Create the random subset by taking the first k elements from the sorted vector
    FullEdges subset;
    for (NodeID i = 0; i < k; i++) {
        subset.push_back(vec[indices[i]]);
    }

    return subset;
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

void longest_path(NodeID n, FullEdges &edges) {
    Graph graph(n, edges);

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

    int longest_path = 0;
    for (int i = 0; i < n; i++)
        longest_path = std::max(longest_path, count[i]);

    std::clog << "Longest Path = " << longest_path << std::endl;
}

GraphInfo create_complete_graph(const ArgsType& args) {
    NodeID n = convertStringToInt<NodeID>(args[0]);
    FullEdges edges;

    for (NodeID i = 0; i < n; i++) {
        for (NodeID j = 0; j < n; j++) {
            if (i != j) {
                edges.emplace_back(i, j, 1);
            }
        }
    }
    return {n, edges};
}

GraphInfo create_cycle_graph(const ArgsType& args) {
    NodeID n = convertStringToInt<NodeID>(args[0]);

    FullEdges edges;

    for (NodeID i = 0; i < n; i++)
        edges.emplace_back(i, (i+1) % n, 1);

    return {n, edges};
}

GraphInfo load_graph(const ArgsType& args) {
    std::string filename = args[0];

    std::clog << "Loading graph from file: " << filename << std::endl;
    Graph graph = readGraph(filename);

    std::clog << "Constructing edges" << std::endl;
    FullEdges edges;
    NodeID n = graph.numberOfNodes();
    for (NodeID v = 0; v < n; v++) {
        for (auto const& e : graph.getEdgesOf(v)) {
            edges.emplace_back(v, e.target, e.weight);
        }
    }
    return {n, edges};
}

// GraphInfo create_random_graph(const ArgsType& args) {
//     NodeID n = convertStringToInt<NodeID>(args[0]);
//     double p = std::stod(args[1]);
//     RandIntGen<int> gen(120);
//
//
//     FullEdges edges;
//     for (NodeID i = 0; i < n; i++) {
//         for (NodeID j = 0; j < n; j++) {
//             if (i != j) {
//                 edges.emplace_back(i, j, gen() - 60);
//             }
//         }
//     }
//
//     // Shuffle the indices randomly
//     std::random_device rd;
//     std::mt19937 g(rd());
//     std::shuffle(edges.begin(), edges.end(), g);
//     edges.resize(static_cast<EdgeID>(p * n * (n - 1)));
//
//     return {n, edges};
// }

GraphInfo create_random_restricted_graph4(const ArgsType& args) {
    NodeID n = convertStringToInt<NodeID>(args[0]);
    double p = std::stod(args[1]);

    // Creating a random graph with all weights = 2
    double num_edges = p * n * (n - 1);
    FullEdges edges = augment_graph(n, {}, p, 2);

    // Maximising number of edges with weight -1
    std::unordered_set<NodeID> to_be_visited;
    for (NodeID i = 0; i < n; i++) {
        to_be_visited.insert(i);
    }

    Distances potentials(n, -3);
    std::vector<std::vector<EdgeID>> node_to_neighbors(n);  // Maps a node to its edges

    for (EdgeID i = 0; i < edges.size(); i++) {
        node_to_neighbors[std::get<0>(edges[i])].push_back(i);
    }

    while (!to_be_visited.empty()) {
        NodeID tmp_source = *std::next(to_be_visited.begin(), std::rand() % to_be_visited.size());

        std::unordered_set<EdgeID> visited_edge_idx;
        std::unordered_set<NodeID> visited_nodes;

        potentials[tmp_source] = 0;
        std::queue<NodeID> queue;
        queue.push(tmp_source);

        while (!queue.empty()) {
            NodeID curr_source = queue.front();
            queue.pop();

            for (EdgeID edge_idx : node_to_neighbors[curr_source]) {
                FullEdge& edge = edges[edge_idx];

                if (to_be_visited.count(std::get<1>(edge)) == 0) {
                    std::get<2>(edge) = -1;
                    continue;
                }

                if (potentials[std::get<1>(edge)] == -3) {
                    potentials[std::get<1>(edge)] = potentials[curr_source] + std::get<2>(edge);
                    queue.push(std::get<1>(edge));
                }

                visited_edge_idx.insert(edge_idx);
            }

            visited_nodes.insert(curr_source);
        }

        for (EdgeID edge_idx : visited_edge_idx) {
            FullEdge& edge = edges[edge_idx];
            std::get<2>(edge) += potentials[std::get<0>(edge)] - potentials[std::get<1>(edge)] - 1;
        }

        for (NodeID node : visited_nodes) {
            to_be_visited.erase(node);
        }
    }

    return {n, edges};
}

GraphInfo create_random_restricted_graph3(const ArgsType& args) {
    NodeID n = convertStringToInt<NodeID>(args[0]);
    double p = std::stod(args[1]);

    ArgsType args2 = {std::to_string(n - 1), args[1]};

    GraphInfo graph = create_random_restricted_graph4(args2);

    std::get<0>(graph) = n;
    FullEdges& edges = std::get<1>(graph);

    for (FullEdge& e : edges) {
        std::get<0>(e) += 1;
        std::get<1>(e) += 1;
    }

    // Adding supersource
    for (NodeID i = 1; i < n; i++) {
        edges.push_back({0, i, 0});
    }

    return graph;
}

void increase_edges(NodeID n, FullEdges &edges) {
    int INCREASE = 400000; // Double the edge values
    RandIntGen<Distance> gen(INCREASE);

    for (auto &[tail, head, weight] : edges)
        if (weight >= 0)
            weight += gen();
}

void increase_edges_coefficient(NodeID n, FullEdges &edges) {
    int INCREASE = 5;
    long max_weight = 0;
    RandIntGen<Distance> gen(INCREASE);
    for (auto &[tail, head, weight] : edges)
        if (weight >= 0) {

            weight += ( gen() * std::abs(static_cast<long long>(tail) - static_cast<long long>(head)));

            max_weight = std::max(max_weight, weight);
        }

    std::clog << "Max Weight = " << max_weight << std::endl;
    longest_path(n, edges);
}

void decrease_edges(NodeID n, FullEdges &edges) {
    int DECREASE = 10000;  // Small decrease
    RandIntGen<Distance> gen(DECREASE);

    for (auto &[tail, head, weight] : edges)
        if (weight < 0)
            weight -= gen();
}

void potential_transformation(NodeID n, FullEdges &edges) {
    long long int MIN_POTENTIAL = -4 * 1664730, MAX_POTENTIAL = 4 * 1664730;
    RandIntGen<Distance> gen(MAX_POTENTIAL - MIN_POTENTIAL);

    Distances potential(n, 0);
    for (NodeID i = 0; i < n; i++)
        potential[i] = MIN_POTENTIAL + gen();

    for (auto &[tail, head, weight] : edges)
        weight += potential[tail] - potential[head];
}

void add_cycle_to_graph(NodeID n, FullEdges &edges, int add_cycle) {
    if (add_cycle == 0)
        return;

    int length = 0;
    if (add_cycle == 1)
        length = 3;

    if (add_cycle == 2)
        length = sqrt(n);

    if (add_cycle == 3)
        length = n - 5;

    std::vector<NodeID> perm(n);
    for (NodeID i = 0; i < n; ++i) perm[i] = i;

    for (NodeID i = n - 1; i > 0; --i) {
        RandIntGen<NodeID> gen(i);
        NodeID j = gen();
        std::swap(perm[i], perm[j]);
    }

    for (NodeID i = 0; i < length; ++i)
        edges.emplace_back(perm[i], perm[(i + 1) % length], static_cast<Distance>(0));

    // edges.emplace_back(perm[length-1], perm[0], static_cast<Distance>(-1));
}

GraphInfo create_random_graph(const ArgsType& args) {
    NodeID n = convertStringToInt<NodeID>(args[0]);
    double p = std::stod(args[1]);

    FullEdges edges = augment_graph(n, {}, p, 0);

    increase_edges_coefficient(n, edges);
    //increase_edges(n, edges);
    potential_transformation(n, edges);
    return {n, edges};
}

GraphInfo create_random_graph_grid(const ArgsType& args) {
    if (args.size() < 2) {
        std::cerr << "Invalid number of arguments for random_graph_grid (need 2: rows cols)" << std::endl;
        return {0, {}};
    }

    NodeID rows = convertStringToInt<NodeID>(args[0]);
    NodeID cols = convertStringToInt<NodeID>(args[1]);

    if (rows == 0 || cols == 0) {
        return {0, {}};
    }

    NodeID n = rows * cols;
    FullEdges edges;
    edges.reserve(4 * n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<Distance> noise_dist(0, 3);
    std::uniform_int_distribution<Distance> slack_dist(0, 4);

    // A smooth downhill potential plus small noise creates many negative edges
    // while keeping every directed cycle positive after adding a per-edge slack.
    constexpr Distance row_scale = 4;
    constexpr Distance col_scale = 4;

    Distances potential(n, 0);
    for (NodeID r = 0; r < rows; ++r) {
        for (NodeID c = 0; c < cols; ++c) {
            NodeID id = r * cols + c;
            potential[id] = static_cast<Distance>((rows - 1 - r) * row_scale + (cols - 1 - c) * col_scale + noise_dist(gen));
        }
    }

    auto add_edge = [&](NodeID u, NodeID v) {
        Distance weight = static_cast<Distance>(1 + slack_dist(gen) + potential[v] - potential[u]);
        edges.emplace_back(u, v, weight);
    };

    for (NodeID r = 0; r < rows; ++r) {
        for (NodeID c = 0; c < cols; ++c) {
            NodeID u = r * cols + c;

            if (r + 1 < rows) {
                NodeID v = (r + 1) * cols + c;
                add_edge(u, v);
                add_edge(v, u);
            }

            if (c + 1 < cols) {
                NodeID v = r * cols + (c + 1);
                add_edge(u, v);
                add_edge(v, u);
            }
        }
    }

    return {n, edges};
}

GraphInfo bad_for_gor(NodeID n) {
    FullEdges edges;

    for (int i = 0; i < n - 1; i++)
        edges.emplace_back(i, i+1, 0);

    RandIntGen<NodeID> gen(n - 1);
    RandIntGen<NodeID> gen2(20);
    for (int e = 0 ; e < n * 10; e++) {
        NodeID i = gen();
        NodeID j = gen();
        if (i < j)
            edges.emplace_back(i, j, j - i + 1 +gen());
        else
            edges.emplace_back(i, j, gen2() * 100);
    }

    return {n, edges};
}

GraphInfo scaling_step(const ArgsType& args) {
    GraphInfo graph = load_graph(args);

    NodeID n = std::get<0>(graph);
    FullEdges& edges = std::get<1>(graph);

    Distance most_negative = 0;

    for (FullEdge& e : edges) {
        if (std::get<2>(e) < most_negative) {
            most_negative = std::get<2>(e);
        }
    }

    most_negative = -most_negative;

    for (FullEdge& e : edges) {
        Distance &weight = std::get<2>(e);
        //  \ceil{3 w(e)/ (W+1)} + 1
        weight = static_cast<Distance>(std::ceil(static_cast<double>(3 * weight) / (most_negative + 1)) + 1);
    }

    return {n, edges};
}

// Count the number of negative edges
EdgeID countNegativeWeightEdges(const FullEdges &edges) {
    EdgeID negative_edges_count = 0;
    for (const auto &edge: edges) {
        if (std::get<2>(edge) < 0) {
            negative_edges_count++;
        }
    }
    return negative_edges_count;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Invalid number of arguments: ./create_graph <graph type>" << std::endl;
        return 1;
    }

    std::string graph_type = argv[1];

    // Handle arguments
    ArgsType args;
    bool mute = false;
    bool do_perm = false;

    // Do augmentation, fraction of edges, fixed weight
    // Only for BAD instances
    std::tuple<bool, double, Distance> augment = {false, 0.0, 0};

    bool e_increase = false;
    bool e_decrease = false;
    bool p_transformation = false;
    int add_cycle = 0;


    for (int i = 2; i < argc; i++) {
        // permutation
        if (std::string(argv[i]) == "-p") {
            do_perm = true;
            continue;
        }
        //augmentation
        if (std::string(argv[i]) == "-a") {
            if (i + 2 >= argc) {
                std::cerr << "Invalid number of arguments for -a (need 2: fraction of edges (less than 1) or multiplicative factor (greater than 1) and fixed weight)" << std::endl;
                return 1;
            }
            augment = {true, std::stod(argv[i + 1]), convertStringToInt<Distance>(argv[i + 2])};
            i += 2;  // the for loop will increase i by an additional 1
            continue;
        }
        //increase edges >= 0
        if (std::string(argv[i]) == "-e") {
            e_increase = true;
            continue;
        }
        //decrease edges < 0
        if (std::string(argv[i]) == "-d") {
            e_decrease = true;
            continue;
        }
        // potential transformation
        if (std::string(argv[i]) == "-t") {
            p_transformation = true;
            continue;
        }
        // add cycle
        if (std::string(argv[i]) == "-c") {
            if (i + 1 >= argc) {
                std::cerr << "Invalid number of arguments for -c (need 1: type of cycle to add)" << std::endl;
                return 1;
            }
            add_cycle = convertStringToInt<int>(argv[i + 1]);
            if (! (add_cycle >= 0 && add_cycle <= 3)) {
                std::cerr << "Invalid cycle adding" << std::endl;
                return 1;
            }
            i++;
            continue;
        }
        // mute
        if (std::string(argv[i]) == "-m") {
            mute = true;
        }
        args.push_back(argv[i]);
    }

    // Map graph type to function
    std::unordered_map<std::string, GraphCreatorType> graph_type_to_function = {
        {"complete_unit_graph", create_complete_graph},
        {"cycle", create_cycle_graph},
        {"load_graph", load_graph},
        {"random_graph_grid", create_random_graph_grid},
        {"random_graph", create_random_graph},
        {"random_restricted_graph3", create_random_restricted_graph3},
        {"random_restricted_graph4", create_random_restricted_graph4},
        {"scaling_step", scaling_step},
    };

    if (graph_type_to_function.count(graph_type) == 0) {
        std::cout << "Invalid graph type: " << graph_type << ". Valid types: " << std::endl;
        for (const auto& [k, v] : graph_type_to_function) {
            std::cout << "\t" << k << std::endl;
        }
        return 1;
    }

    // Additional operations
    GraphInfo nodes_edges = graph_type_to_function[graph_type](args);
    auto& [nodes, edges] = nodes_edges;

    if (std::get<0>(augment)) {
        // std::cout << "Augmenting " << std::get<1>(augment) << " " << std::get<2>(augment) << std::endl << std::flush;
        auto new_edges = augment_graph(nodes, edges, std::get<1>(augment), std::get<2>(augment));
        edges.insert(edges.end(), new_edges.begin(), new_edges.end());
    }

    if (e_increase) increase_edges(nodes, edges);
    if (e_decrease) decrease_edges(nodes, edges);

    add_cycle_to_graph(nodes, edges, add_cycle);

    if (p_transformation) potential_transformation(nodes, edges);



    if (do_perm) {
        std::vector<NodeID> node_map(nodes);
        std::iota(node_map.begin(), node_map.end(), 0);
        shuffle_vec(node_map.begin(), node_map.end());
        for (auto & [tail, head, weight] : edges) {
            tail = node_map[tail];
            head = node_map[head];
        }
    }

    std::sort(edges.begin(), edges.end());
    Graph graph(nodes, edges);

    std::clog << "Done." << std::endl << std::flush;
    std::clog << "Number of negative weight edges: " << countNegativeWeightEdges(edges) << std::endl << std::flush;
    std::clog << "Number of SCC: " << decomposeIntoSCCs(graph).size() << std::endl << std::flush;

    if (!mute) graph.format_print();
    return 0;
}