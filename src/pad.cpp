//
// Created by adumi on 1/27/26.
//

#include "pad.h"

#include "measurement_tool.h"

Graph get_new_graph(const Graph &graph, const std::vector<bool> &u) {
    NodeID n = graph.numberOfNodes();
    std::vector<NodeID> original_id, new_id(n);
    for (int i = 0; i < n; i++)
        if (u[i]) {
            new_id[i] = original_id.size();
            original_id.push_back(i);
        }

    std::vector<FullEdge> new_edges;

    for (int i = 0; i < n; i++)
        if (u[i]) {
            for (const auto edge: graph.getEdgesOf(i))
                if (u[edge.target])
                    new_edges.emplace_back(new_id[i], new_id[edge.target], edge.weight);
        }

    return Graph(original_id.size(), new_edges, original_id);
}

EdgeID grow_ball(const Distance new_d, int &start, const int end, const std::vector<Distance> &d,
                 const std::vector<NodeID> &order,
                 const std::vector<EdgeID> &vol, const Distance diameter) {

    if (end == 0) // NO ELEMENTS IN THE BALL YET
        return 0;

    const Distance delta = diameter / config::pad_rounds;

    if (new_d == d[order[end]])
        return 0;

    while (d[order[start]] + delta < new_d) {
        if (start == end) // this means that no padding
            return vol[start];

        if (d[order[start]] != d[order[start + 1]]) {
            // included everything in ball
            if ((config::pad_alpha + 1) * vol[start] >= config::pad_alpha * vol[end])
                return vol[start];
        }

        start++;
    }

    return 0;
}

std::pair<int, int> grow_ball_heavy(const Graph &graph, Distance diameter, std::vector<NodeID> &order, NodeID s,
                                    Orientation orientation) {
    NodeID n = graph.numberOfNodes();
    int start = 0, end = 0;
    pad::GraphHeap q(n);
    Distances d(n, c::infty);
    std::vector<EdgeID> volumes(n + 1, 0);

    d[s] = 0;

    q.insert(s, 0);

    while (!q.empty()) {
        if (q.minKey() > diameter / 10)
            break;

        NodeID from;
        Distance dist;
        q.deleteMin(from, dist);

        start++;
        end++;
        order[end] = from;
        volumes[end] = volumes[end - 1] + graph.getAllDegreeOf(from);

        // Dijkstra normal graph
        for (auto const &edge: graph.getEdgesOf(from, orientation)) {
            auto tentative_dist = d[from] + std::max(static_cast<Distance>(0), edge.weight);
            if (tentative_dist < d[edge.target]) {
                d[edge.target] = tentative_dist;
                if (q.contains(edge.target)) {
                    q.decreaseKey(edge.target, tentative_dist);
                } else {
                    q.insert(edge.target, tentative_dist);
                }
            }
        }
    }

    EdgeID ball = 0;
    while (!q.empty()) {
        NodeID from;
        Distance dist;
        q.deleteMin(from, dist);

        ball = grow_ball(dist, start, end, d, order, volumes, diameter);

        if (ball)
            break;

        end++;
        order[end] = from;
        volumes[end] = volumes[end - 1] + graph.getAllDegreeOf(from);

        // Dijkstra normal graph
        for (auto const &edge: graph.getEdgesOf(from, orientation)) {
            auto tentative_dist = d[from] + std::max(static_cast<Distance>(0), edge.weight);
            if (tentative_dist < d[edge.target]) {
                d[edge.target] = tentative_dist;
                if (q.contains(edge.target)) {
                    q.decreaseKey(edge.target, tentative_dist);
                } else {
                    q.insert(edge.target, tentative_dist);
                }
            }
        }
    }

    if (!ball)
        grow_ball(c::infty, start, end, d, order, volumes, diameter);


    //TODO() remove check //////////////////////////////////////////////////////////////////////////////////////////////
    if (volumes[start] < graph.numberOfEdges() || d[order[end]] > diameter / 5 + 1) {
        PRINT("INCORRECT HEAVY LIMITS");
        exit(-1);
    }
    std::vector<bool> u(n, false), pad(n, false);
    for (int i = 1; i <= start; i++)
        u[order[i]] = pad[order[i]] = true;
    for (int i = start + 1; i <= end; i++)
        pad[order[i]] = true;
    padding_check(graph, u, pad, diameter, orientation);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    return {start, end};
}

std::vector<Graph> padded_decomposition_heavy(const Graph &graph, Distance diameter, NodeID s) {
    // under diameter/10 I found a ball of big volume
    PRINT("HEAVY");
    NodeID n = graph.numberOfNodes();

    std::vector<NodeID> order_plus(n + 1, 0), order_minus(n + 1, 0);

    auto [start_plus, end_plus] = grow_ball_heavy(graph, diameter, order_plus, s, Orientation::OUT);
    auto [start_minus, end_minus] = grow_ball_heavy(graph, diameter, order_minus, s, Orientation::IN);

    std::vector<Graph> X(3);

    std::vector<int> intersect(n);
    std::vector<bool> u(n, false);

    for (int i = 1; i <= end_plus; i++)
        intersect[order_plus[i]]++;
    for (int i = 1; i <= end_minus; i++)
        intersect[order_minus[i]]++;
    for (int i = 0; i < n; i++)
        if (intersect[i] == 2)
            u[i] = true;
    X[2] = get_new_graph(graph, u);

    u.assign(u.size(), false);
    for (int i = 1; i <= end_plus; i++)
        u[order_plus[i]] = true;
    for (int i = 1; i <= start_minus; i++)
        u[order_minus[i]] = false;
    X[0] = get_new_graph(graph, u);

    u.assign(u.size(), true);
    for (int i = 1; i <= start_plus; i++)
        u[order_plus[i]] = false;
    X[1] = get_new_graph(graph, u);

    return X;
}

std::vector<Graph> padded_decomposition(Graph &graph, Distance diameter) {
    NodeID n = graph.numberOfNodes();
    EdgeID m = graph.numberOfEdges();

    std::vector<bool> u_plus(n, false);
    std::vector<bool> u_minus(n, false);
    std::vector<bool> pad_plus(n, false);
    std::vector<bool> pad_minus(n, false);

    EdgeID vol_u_plus = 0;
    EdgeID vol_u_minus = 0;

    pad::GraphHeap q_plus(n);
    pad::GraphHeap q_minus(n);

    Distances d_plus(n, c::infty);
    Distances d_minus(n, c::infty);

    std::vector<NodeID> order_plus(n + 1, 0);
    std::vector<NodeID> order_minus(n + 1, 0);

    std::vector<EdgeID> volumes_plus(n + 1, 0);
    std::vector<EdgeID> volumes_minus(n + 1, 0);

    NodeID s = 0;
    while (vol_u_plus * 2 < m && vol_u_minus * 2 < m) {
        // Finding initial node
        if (u_plus[s] || u_minus[s]) {
            s++;
            continue;
        }

        d_plus[s] = 0;
        d_minus[s] = 0;

        q_plus.insert(s, 0);
        q_minus.insert(s, 0);

        int start_plus = 1, start_minus = 1, end = 0;
        EdgeID ball_plus = 0, ball_minus = 0;

        while (true) {
            // I always find the ball...
            bool plus_not_done = false, minus_not_done = false;

            if (!q_plus.empty() && ball_plus == 0) {
                plus_not_done = true;
                ball_plus = grow_ball(q_plus.minKey(), start_plus, end, d_plus, order_plus, volumes_plus, diameter);
            } else {
                if (ball_plus == 0)
                    ball_plus = grow_ball(c::infty, start_plus, end, d_plus, order_plus, volumes_plus, diameter);
            }


            if (!q_minus.empty() && ball_minus == 0) {
                minus_not_done = true;
                ball_minus = grow_ball(q_minus.minKey(), start_minus, end, d_minus, order_minus, volumes_minus, diameter);
            } else {
                if (ball_minus == 0)
                    ball_minus = grow_ball(c::infty, start_minus, end, d_minus, order_minus, volumes_minus, diameter);
            }

            //TODO remove check ////////////////////////////////////////////////////////////////////////////////////////
            if (ball_plus && d_plus[order_plus[start_plus]] > diameter / 10) {
                PRINT("WRONG SMALL LIMIT");
                exit(-1);
            }
            if (ball_minus && d_minus[order_minus[start_minus]] > diameter / 10) {
                PRINT("WRONG SMALL LIMIT");
                exit(-1);
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (ball_plus && ball_minus)
                break;

            if ((ball_plus || ball_minus) && (ball_plus < m && ball_minus < m))
                break;


            end++;

            if (plus_not_done) {
                NodeID from_plus; Distance dist_plus;
                q_plus.deleteMin(from_plus, dist_plus);
                order_plus[end] = from_plus;
                volumes_plus[end] = volumes_plus[end - 1] + graph.getAllDegreeOf(from_plus);

                // Dijkstra normal graph
                for (auto const &edge: graph.getEdgesOf(from_plus)) {
                    if (u_plus[edge.target])
                        continue;

                    auto tentative_dist = d_plus[from_plus] + std::max(static_cast<Distance>(0), edge.weight);
                    if (tentative_dist < d_plus[edge.target]) {
                        d_plus[edge.target] = tentative_dist;

                        if (q_plus.contains(edge.target)) {
                            q_plus.decreaseKey(edge.target, tentative_dist);
                        } else {
                            q_plus.insert(edge.target, tentative_dist);
                        }
                    }
                }
            }

            if (minus_not_done) {
                NodeID from_minus; Distance dist_minus;
                q_minus.deleteMin(from_minus, dist_minus);
                order_minus[end] = from_minus;
                volumes_minus[end] = volumes_minus[end - 1] + graph.getAllDegreeOf(from_minus);

                // For reversed graph
                for (auto const &edge: graph.getEdgesOf(from_minus, Orientation::IN)) {
                    if (u_minus[edge.target])
                        continue;

                    auto tentative_dist = d_minus[from_minus] + std::max(static_cast<Distance>(0), edge.weight);
                    if (tentative_dist < d_minus[edge.target]) {
                        d_minus[edge.target] = tentative_dist;

                        if (q_minus.contains(edge.target)) {
                            q_minus.decreaseKey(edge.target, tentative_dist);
                        } else {
                            q_minus.insert(edge.target, tentative_dist);
                        }
                    }
                }
            }
        }

        // I should definitely not be here never ever...
        if (ball_plus == 0 && ball_minus == 0) {
            exit(-1);
        }

        if (ball_plus > m && ball_minus > m)
            return padded_decomposition_heavy(graph, diameter, s);

        if (ball_plus == 0) ball_plus = 2 * m + 1;
        if (ball_minus == 0) ball_minus = 2 * m + 1;

        if (ball_plus <= ball_minus) {
            vol_u_plus += ball_plus;
            // update u_plus and padding vectors
            for (int i = 1; i <= start_plus; i++) {
                u_plus[order_plus[i]] = true;
                pad_plus[order_plus[i]] = true;
            }
            for (int i = start_plus + 1; i <= end; i++)
                if (d_plus[order_plus[start_plus]] + diameter / config::pad_rounds >= d_plus[order_plus[i]])
                    pad_plus[order_plus[i]] = true;
        } else {
            vol_u_minus += ball_minus;
            // update u_minus and padding vectors
            for (int i = 1; i <= start_minus; i++) {
                u_minus[order_minus[i]] = true;
                pad_minus[order_minus[i]] = true;
            }
            for (int i = start_minus + 1; i <= end; i++)
                if (d_minus[order_minus[start_minus]] + diameter / config::pad_rounds >= d_minus[order_minus[i]])
                    pad_minus[order_minus[i]] = true;
        }


        // reset all
        for (int i = 1; i <= end; i++) {
            d_plus[order_plus[i]] = d_minus[order_minus[i]] = c::infty;
            volumes_plus[i] = volumes_minus[i] = 0;
            order_plus[i] = order_minus[i] = 0;
        }
        while (!q_plus.empty()) {
            NodeID v;
            Distance d;
            q_plus.deleteMin(v, d);
            d_plus[v] = c::infty;
        }

        while (!q_minus.empty()) {
            NodeID v;
            Distance d;
            q_minus.deleteMin(v, d);
            d_minus[v] = c::infty;
        }
    }

    // We exited the loop = build the graphs
    std::vector<Graph> X(2);
    std::vector<bool> u, pad;
    if (vol_u_plus > vol_u_minus) {
        u = std::move(u_plus);
        pad = std::move(pad_plus);
        PRINT("LIGHT - PLUS");
        //TODO remove check ////////////////////////////////////////////////////////////////////////////////////////////
        padding_check(graph, u, pad, diameter, Orientation::OUT);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } else {
        u = std::move(u_minus);
        pad = std::move(pad_minus);
        PRINT("LIGHT - MINUS");
        //TODO remove check ////////////////////////////////////////////////////////////////////////////////////////////
        padding_check(graph, u, pad, diameter, Orientation::IN);
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

    for (NodeID i = 0; i < n; i++)
        u[i] = (!u[i]);

    X[0] = get_new_graph(graph, pad);
    X[1] = get_new_graph(graph, u);

    return X;
}

std::optional<Distances> pad::PADAlg::runMainAlg(Graph &graph, Distance diameter, int level) {
    const NodeID n = graph.numberOfNodes();
    const EdgeID m = graph.numberOfEdges();

    stats.max_recursion_level = std::max(level, stats.max_recursion_level);

    if (n <= config::pad_small) {
        Distances p(n, 0);
        stats.in_padding = false;
        if (config::cycle_detection)
            return runLazyDijkstra(graph, p, diameter, n + 2);

        return runLazyDijkstra(graph, p, diameter, -1); // NO CYCLE DETECTION = GO UNTIL END
    }

    // Compute 2-approx of the diameter - ONE TIME I USE BCF
    const Distances &d_out = bcf::runDijkstra(graph, 0, c::infty, Orientation::OUT);
    const Distances &d_in = bcf::runDijkstra(graph, 0, c::infty, Orientation::IN);

    diameter = std::min(
        diameter, *std::max_element(d_out.begin(), d_out.end()) + *std::max_element(d_in.begin(), d_in.end()));

    PRINT("Recursion Level: " << level);
    if (diameter <= config::pad_rounds) {
        PRINT("END OF RECURSION - SMALL DIAMETER: " << diameter);
        Distances p(n, 0);
        if (config::cycle_detection)
            return runLazyDijkstra(graph, p, diameter, diameter + 2);

        return runLazyDijkstra(graph, p, diameter, -1);
    }

    PRINT("PADDED DECOMPOSITION: n = " << n << ", m = " << m << ", diameter = " << diameter);
    auto X = padded_decomposition(graph, diameter);

    NodeID H_n = 0;
    for (const auto &i: X)
        H_n += i.numberOfNodes();

    Distances phi(H_n);
    std::vector<std::vector<NodeID> > membership(n);
    std::vector<NodeID> global_id(H_n);
    std::vector<NodeID> aux_n(X.size(), 0);

    aux_n[1] = X[0].numberOfNodes();
    if (X.size() == 3) {
        stats.heavy_calls++;
        aux_n[2] = aux_n[1] + X[1].numberOfNodes();
    }

    for (int i = 0; i < X.size(); i++) {
        Distances potential(X[i].numberOfNodes());
        auto components = decomposeIntoSCCs(X[i]);
        Distance cDiameter = (i < 2) ? diameter : (diameter / 2 + 1);

        PRINT("RECURSING ON - " << components.size() << " PROBLEMS");
        for (auto &component: components) {
            auto opt_comp_potential = runMainAlg(component, cDiameter, level + 1);

            if (!opt_comp_potential.has_value())
                return {};

            auto component_potential = std::move(opt_comp_potential.value());

            for (NodeID j = 0; j < component.numberOfNodes(); j++)
                potential[component.global_id[j]] = component_potential[j];
        }

        bcf::fixDagEdges(X[i], components, potential);

        for (NodeID j = 0; j < X[i].numberOfNodes(); j++) {
            NodeID H_id = aux_n[i] + j;
            NodeID G_id = X[i].global_id[j];

            phi[H_id] = potential[j];
            membership[G_id].push_back(H_id);
            global_id[H_id] = G_id;
        }
    }

    PRINT("ADDING EDGES TO H");
    bool has_padding = false;
    std::vector<FullEdge> e;
    int negative_edges = 0;
    for (int i = 0; i < n; i++) {
        if (membership[i].size() != 1)
            has_padding = true;

        for (auto edge: graph.getEdgesOf(i))
            for (auto u: membership[i])
                for (auto v: membership[edge.target]) {
                    e.emplace_back(u, v, edge.weight);
                    if (edge.weight + phi[u] - phi[v] < 0)
                        negative_edges++;
                }
    }

    MEASUREMENT::addInt(EXP::NEGATIVE_EDGES_IN_DECOMPOSITION, negative_edges);

    PRINT("CREATING H");
    stats.decomposition_calls++;
    if (has_padding) {
        PRINT("HAS PADDING");
        stats.decomposition_calls_with_padding++;
    } else
        PRINT("NO PADDING");

    auto H = Graph(H_n, e);
    PRINT("DONE CREATING");

    PRINT("    RUNNING LAZY DIJKSTRA: H_n = " << H_n << ", e.size() = " << e.size() << ", n = " << n << ", m = " << m <<
        ", diameter = " << diameter);
    std::optional<Distances> optional_H_potential;
    if (config::pad_use_lazy) {
        // A max number of rounds OR stopping condition from scaling
        stats.in_padding = true;
        if (config::cycle_detection)
            optional_H_potential = runLazyDijkstra(H, phi, diameter, -1);
        else
            optional_H_potential = runLazyDijkstra(H, phi, diameter, -1);
    } else
        optional_H_potential = gor(H, phi);

    if (!optional_H_potential.has_value())
        return {};
    auto H_potential = std::move(optional_H_potential.value());

    Distances potential(n);
    for (int i = 0; i < H_n; i++)
        potential[global_id[i]] = H_potential[i] + phi[i];

    for (NodeID from = 0; from < graph.numberOfNodes(); from++) {
        for (auto const &edge: graph.getEdgesOf(from)) {
            auto pot_edge_w = edge.weight + potential[from] - potential[edge.target];
            if (pot_edge_w < 0) {
                PRINT("SUPER BAD!!!");
                exit(-1);
            }
        }
    }
    PRINT("DONE LAZY DIJKSTRA");
    return potential;
}



std::optional<Distances> pad::runLazyDijkstra(const Graph &graph, const Distances &potential, Distance diameter,
                                              int max_rounds) {
    NodeID n = graph.numberOfNodes();
    Distances distance(n, c::infty);
    Distances positive(n, 0);
    std::vector<NodeID> bellman_phase;
    GraphHeap q(n);
    int rounds = 0;


    bellman_phase.reserve(n);
    // This is done using potentials
    for (NodeID i = 0; i < n; i++) {
        distance[i] = -potential[i];
        q.insert(i, distance[i]);
    }

    while (!q.empty()) {
        // PRINT("ROUNDS: " << rounds);
        if (rounds == max_rounds) {
            PRINT("MAX ROUNDS REACHED: " << max_rounds);
            return {};
        }

        // run Dijkstra phase
        bellman_phase.clear();
        while (!q.empty()) {
            Distance dist;
            NodeID from;
            q.deleteMin(from, dist);

            if (dist > distance[from]) continue;

            if (distance[from] + potential[from] < 0 && positive[from] > (config::pad_scaling_factor - 1) * diameter) {
                if (stats.in_padding)
                    MEASUREMENT::addInt(EXP::LAZY_IN_PADDING, rounds);
                else
                    MEASUREMENT::addInt(EXP::LAZY_IN_SMALL, rounds);

                PRINT("    HEURISTIC WORKS");
                return {};
            }

            bellman_phase.emplace_back(from);

            for (auto const &edge: graph.getEdgesOf(from)) {
                Distance weight = edge.weight + potential[from] - potential[edge.target];

                if (weight < 0)
                    continue;

                auto tentative_dist = distance[from] + weight;
                if (tentative_dist < distance[edge.target]) {
                    distance[edge.target] = tentative_dist;
                    positive[edge.target] = positive[from] + std::max(static_cast<Distance>(0), edge.weight);

                    q.insert(edge.target, tentative_dist);
                }
            }
        }

        // Relax all negative edges from vertices explored in dijkstra phase
        for (auto const &from: bellman_phase)
            for (auto const &edge: graph.getEdgesOf(from)) {
                auto weight = edge.weight + potential[from] - potential[edge.target];

                //if (weight >= 0) continue;

                auto tentative_dist = distance[from] + weight;
                if (tentative_dist <= distance[edge.target]) {
                    if (tentative_dist < distance[edge.target]) {
                        distance[edge.target] = tentative_dist;
                        q.insert(edge.target, tentative_dist);

                        positive[edge.target] = positive[from] + std::max(static_cast<Distance>(0), edge.weight);
                    } else {
                        positive[edge.target] = std::max(positive[edge.target],
                                                         positive[from] + std::max(
                                                             static_cast<Distance>(0), edge.weight));
                    }
                }
            }

        rounds++;
    }

    if (stats.in_padding)
        MEASUREMENT::addInt(EXP::LAZY_IN_PADDING, rounds);
    else
        MEASUREMENT::addInt(EXP::LAZY_IN_SMALL, rounds);

    if (rounds > 1)
        PRINT("    rounds of lazy dijkstra: " << rounds);

    return distance;
}

std::optional<Distances> pad::scaling_early_finish(const Graph &graph, const Graph &current_graph, NodeID source) {
    NodeID n = graph.numberOfNodes();
    // Dijkstra for the tree
    std::vector<NodeID> p(n, -1);
    Distances distances(n, c::infty);
    distances[source] = 0;
    GraphHeap q(n);
    q.insert(source, 0);

    while (!q.empty()) {
        Distance dist;
        NodeID from;
        q.deleteMin(from, dist);
        if (dist > distances[from]) continue;
        for (auto const &edge: current_graph.getEdgesOf(from)) {
            auto tentative_dist = distances[from] + std::max(static_cast<Distance>(0), edge.weight);
            if (tentative_dist < distances[edge.target]) {
                p[edge.target] = from;
                distances[edge.target] = tentative_dist;
                q.insert(edge.target, tentative_dist);
            }
        }
    }

    // Compute the final distances on the original graph
    std::fill(distances.begin(), distances.end(), c::infty);

    std::vector<std::vector<std::pair<NodeID, Distance> > > adj(n);
    for (NodeID v = 0; v < n; v++)
        for (auto e: graph.getEdgesOf(v))
            if (v == p[e.target])
                distances[e.target] = std::min(distances[e.target], e.weight);

    for (NodeID v = 0; v < n; v++)
        if (p[v] != -1)
            adj[p[v]].emplace_back(v, distances[v]);

    // final DFS
    distances[source] = 0;
    std::vector<NodeID> st;
    st.reserve(n);
    st.push_back(source);

    while (!st.empty()) {
        NodeID v = st.back();
        st.pop_back();

        for (const auto &[u, w]: adj[v]) {
            distances[u] = distances[v] + w;
            st.push_back(u);
        }
    }

    if (isResultCorrect(graph, distances, source))
        return distances;

    return {};
}

bool pad::fast_admissible_graph_check(const Graph &graph, const Distances &potential) {
    NodeID n = graph.numberOfNodes();

    int scc_count = 0;
    std::vector<int> scc(n, 0); // use as visitation vector (scc[i] = 0) means node visited
    std::vector<bool> admissible(n, false);

    std::vector<EdgeRange::iterator> Current(n), Current_T(n);
    std::vector<EdgeRange> edgesOf(n), edgesOf_T(n);

    for (NodeID i = 0; i < n; i++) {
        edgesOf[i] = graph.getEdgesOf(i, Orientation::OUT);
        edgesOf_T[i] = graph.getEdgesOf(i, Orientation::IN);

        Current[i] = edgesOf[i].begin();
        Current_T[i] = edgesOf_T[i].begin();
    }

    std::vector<int> topo_sort, dfs;
    dfs.reserve(n);
    topo_sort.reserve(n);

    for (NodeID i = 0; i < n; i++) {
        for (auto &edge: graph.getEdgesOf(i))
            if (potential[i] + edge.weight <= potential[edge.target]) {
                admissible[i] = true;
                break;
            }
    }

    for (NodeID i = 0; i < n; i++)
        if (admissible[i] && scc[i] == 0) {
            // not yet visited...
            dfs.push_back(i);

            while (dfs.size() != 0) {
                NodeID v = dfs.back();
                dfs.pop_back();

                bool found = false;
                for (auto arc = Current[v]; arc != edgesOf[v].end() && !found; ++arc) {
                    NodeID u = arc->target;

                    if (!admissible[u] || scc[u] != 0) // Not admissible or visited
                        continue;

                    // if admissible edge
                    if (potential[v] + arc->weight <= potential[u]) {
                        Current[v] = arc;
                        ++Current[v];

                        dfs.push_back(v);

                        scc[u] = 1;
                        dfs.push_back(u);

                        found = true;
                    }
                }

                if (!found) {
                    // t_out
                    topo_sort.push_back(v);
                }
            }
        }

    // reset scc vector
    for (auto it: topo_sort)
        scc[it] = 0;

    for (auto it: std::views::reverse(topo_sort)) {
        if (scc[it] != 0)
            continue;

        scc_count++;
        scc[it] = scc_count;
        dfs.push_back(it);

        while (dfs.size() != 0) {
            NodeID v = dfs.back();
            dfs.pop_back();

            for (auto arc = Current_T[v]; arc != edgesOf_T[v].end(); ++arc) {
                NodeID u = arc->target;

                if (!admissible[u])
                    continue;

                if (potential[u] + arc->weight <= potential[v]) {
                    bool negative = (potential[u] + arc->weight <= potential[v]);

                    if (scc[u] == scc[v] && negative)
                        return true;

                    if (scc[u] == 0) {
                        if (negative)
                            return true;

                        Current_T[v] = ++arc;
                        dfs.push_back(v);

                        scc[u] = scc_count;
                        dfs.push_back(u);

                        break;
                    }
                }
            }
        }
    }

    std::vector<NodeID> scc_size(scc_count, 0);
    NodeID maxx = 0;
    for (int i = 0; i < n; i++)
        if (scc[i])
            scc_size[scc[i] - 1]++;

    for (int i = 0; i < scc_count; i++)
        maxx = std::max(maxx, scc_size[i]);

    MEASUREMENT::addInt(EXP::SCC_ADMISSIBLE_GRAPH, maxx);
    return false;
}

// pad

void padding_check(const Graph &graph, const std::vector<bool> &u, const std::vector<bool> &pad, Distance diameter, Orientation orientation) {
    NodeID n = graph.numberOfNodes();
    bcf::GraphHeap q(n);
    // Simple checks:
    EdgeID vU = 0, vP = 0;
    for (int i = 0; i < n; i++) {
        if (u[i] && !pad[i]) {
            PRINT("NODE IN CORE BUT NOT IN PADDING: BUG1");
            exit(-1);
        }

        if (u[i])
            vU += graph.getAllDegreeOf(i);
        if (pad[i])
            vP += graph.getAllDegreeOf(i);
    }

    if (!((config::pad_alpha * vP) <= (config::pad_alpha + 1) * vU)) {
        PRINT("PADDING RATIO NOT OK: BUG2");
        exit(-1);
    }

    Distances d(n, c::infty);
    for (int i = 0; i < n; i++)
        if (u[i]) {
            d[i] = 0;
            q.insert(i, 0);
        }

    while (!q.empty()) {
        NodeID from; Distance dist;
        q.deleteMin(from, dist);

        if (d[from] < dist)
            continue;

        for (auto edge: graph.getEdgesOf(from, orientation)) {
            Distance weight = std::max(static_cast<Distance>(0), edge.weight);

            Distance tentative = d[from] + weight;
            if (tentative < d[edge.target]) {
                d[edge.target] = tentative;
                q.insert(edge.target, d[edge.target]);
            }
        }
    }

    for (int i = 0; i < n; i++)
        if (!pad[i] && d[i] < diameter / config::pad_rounds) {
            PRINT("PADDING DISTANCE NOT OK: BUG3");
            exit(-1);
        }
}

void check_cycle_correctness() {
    //TODO
}