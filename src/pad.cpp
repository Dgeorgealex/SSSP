//
// Created by adumi on 1/27/26.
//

#include "pad.h"

#include "bcf.h"

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
    pad::GraphHeap q(n + 1);
    Distances d(n, c::infty);
    std::vector<EdgeID> volumes(n, 0);

    d[s] = 0;
    order[0] = s;
    volumes[0] = graph.getDegreeOf(s);

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
        volumes[end] = volumes[end - 1] + graph.getDegreeOf(from);

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
        volumes[end] = volumes[end - 1] + graph.getDegreeOf(from);

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

    return {start, end};
}

std::vector<Graph> padded_decomposition_heavy(const Graph &graph, Distance diameter, NodeID s) {
    // under d/10 I found a ball of big volume
    PRINT("HEAVY");
    NodeID n = graph.numberOfNodes();

    std::vector<NodeID> order_plus(n), order_minus(n);

    auto [start_plus, end_plus] = grow_ball_heavy(graph, diameter, order_plus, s, Orientation::OUT);
    auto [start_minus, end_minus] = grow_ball_heavy(graph, diameter, order_minus, s, Orientation::IN);

    std::vector<Graph> X(3);

    std::vector<int> intersect(n);
    std::vector<bool> u(n, false);

    for (int i = 0; i <= end_plus; i++)
        intersect[order_plus[i]]++;
    for (int i = 0; i <= end_minus; i++)
        intersect[order_minus[i]]++;
    for (int i = 0; i < n; i++)
        if (intersect[i] == 2)
            u[i] = true;
    X[2] = get_new_graph(graph, u);

    u.assign(u.size(), false);
    for (int i = 0; i <= end_plus; i++)
        u[order_plus[i]] = true;
    for (int i = 0; i <= start_minus; i++)
        u[order_minus[i]] = false;
    X[0] = get_new_graph(graph, u);

    u.assign(u.size(), true);
    for (int i = 0; i <= start_plus; i++)
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

    std::vector<NodeID> order_plus(n, 0);
    std::vector<NodeID> order_minus(n, 0);

    std::vector<EdgeID> volumes_plus(n, 0);
    std::vector<EdgeID> volumes_minus(n, 0);

    int start_plus = 0;
    int start_minus = 0;
    int end = 0;

    NodeID s = 0;
    while (vol_u_plus * 2 < m && vol_u_minus * 2 < m) {
        // Finding initial node
        if (u_plus[s] || u_minus[s]) {
            s++;
            continue;
        }

        d_plus[s] = 0;
        d_minus[s] = 0;

        order_plus[0] = s;
        order_minus[0] = s;

        volumes_plus[0] = graph.getDegreeOf(s);
        volumes_minus[0] = graph.getDegreeOf(s);

        q_plus.insert(s, 0);
        q_minus.insert(s, 0);

        EdgeID ball_plus = 0, ball_minus = 0;

        while (true) {
            // I always find the ball...
            NodeID from_plus, from_minus;
            bool plus_not_done = false, minus_not_done = false;


            if (!q_plus.empty()) {
                Distance dist_plus;
                plus_not_done = true;
                q_plus.deleteMin(from_plus, dist_plus);
                ball_plus = grow_ball(dist_plus, start_plus, end, d_plus, order_plus, volumes_plus, diameter);
            } else {
                if (ball_plus == 0)
                    ball_plus = grow_ball(c::infty, start_plus, end, d_plus, order_plus, volumes_plus, diameter);
            }


            if (!q_minus.empty()) {
                Distance dist_minus;
                minus_not_done = true;
                q_minus.deleteMin(from_minus, dist_minus);
                ball_minus = grow_ball(dist_minus, start_minus, end, d_minus, order_minus, volumes_minus, diameter);
            } else {
                if (ball_minus == 0)
                    ball_minus = grow_ball(c::infty, start_minus, end, d_minus, order_minus, volumes_minus, diameter);
            }


            if (ball_plus && ball_minus)
                break;

            if ((ball_plus || ball_minus) && (ball_plus < m && ball_minus < m))
                break;


            end++;


            if (plus_not_done) {
                order_plus[end] = from_plus;
                volumes_plus[end] = volumes_plus[end - 1] + graph.getDegreeOf(from_plus);

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
                order_minus[end] = from_minus;
                volumes_minus[end] = volumes_minus[end - 1] + graph.getDegreeOf(from_minus);

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

        if (ball_plus >= m && ball_minus >= m)
            return padded_decomposition_heavy(graph, diameter, s);

        if (ball_plus == 0) ball_plus = m;
        if (ball_minus == 0) ball_minus = m;

        if (ball_plus <= ball_minus) {
            vol_u_plus += ball_plus;
            // update u_plus and padding vectors
            for (int i = 0; i <= start_plus; i++) {
                u_plus[order_plus[i]] = true;
                pad_plus[order_plus[i]] = true;
            }
            for (int i = start_plus + 1; i <= end; i++)
                if (d_plus[order_plus[start_plus]] + diameter / config::pad_rounds >= d_plus[order_plus[i]])
                    pad_plus[order_plus[i]] = true;
        } else {
            vol_u_minus += ball_minus;
            // update u_minus and padding vectors
            for (int i = 0; i <= start_minus; i++) {
                u_minus[order_minus[i]] = true;
                pad_minus[order_minus[i]] = true;
            }
            for (int i = start_minus + 1; i <= end; i++)
                if (d_minus[order_minus[start_minus]] + diameter / config::pad_rounds >= d_minus[order_minus[i]])
                    pad_minus[order_minus[i]] = true;
        }


        // reset all
        for (int i = 0; i <= end; i++) {
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

        start_plus = start_minus = end = 0;
    }

    // We exited the loop = build the graphs
    PRINT("LIGHT");
    std::vector<Graph> X(2);
    std::vector<bool> u, pad;
    if (vol_u_plus > vol_u_minus) {
        u = std::move(u_plus);
        pad = std::move(pad_plus);
    } else {
        u = std::move(u_minus);
        pad = std::move(pad_minus);
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

    if (n <= config::pad_small) {
        Distances p(n, 0);
        return bcf::runLazyDijkstra(graph, p);
    }

    // Compute 2-approx of the diameter
    const Distances &d_out = bcf::runDijkstra(graph, 0, c::infty, Orientation::OUT);
    const Distances &d_in = bcf::runDijkstra(graph, 0, c::infty, Orientation::IN);

    diameter = std::min(
        diameter, *std::max_element(d_out.begin(), d_out.end()) + *std::max_element(d_in.begin(), d_in.end()));

    PRINT("Recursion Level: " << level);
    if (diameter <= config::pad_rounds) {
        PRINT("END OF RECURSION - SMALL DIAMETER: " << diameter);
        Distances p(n, 0);
        return bcf::runLazyDijkstra(graph, p, Orientation::OUT, diameter + 2);
    }

    PRINT("PADDED DECOMPOSITION: n = " << n << ", m = " << m << ", diameter = " << diameter);
    auto X = padded_decomposition(graph, diameter);

    NodeID H_n = 0;
    for (const auto & i : X)
        H_n += i.numberOfNodes();

    Distances phi(H_n);
    std::vector<std::vector<NodeID> > membership(n);
    std::vector<NodeID> global_id(H_n);
    std::vector<NodeID> aux_n(X.size(), 0);

    aux_n[1] = X[0].numberOfNodes();
    for (int i = 0; i < 2; i++) {
        Distances potential(X[i].numberOfNodes());
        auto components = decomposeIntoSCCs(X[i]);

        PRINT("RECURSING ON - " << components.size() << " PROBLEMS");
        for (auto &component: components) {
            auto opt_comp_potential = runMainAlg(component, diameter, level + 1);

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

    if (X.size() == 3) {
        auto opt_potential = runMainAlg(X[2], diameter / 2, level + 1);

        if (!opt_potential.has_value())
            return {};

        auto potential = std::move(opt_potential.value());

        aux_n[2] = X[0].numberOfNodes() + X[1].numberOfNodes();

        for (NodeID j = 0; j < X[2].numberOfNodes(); j++) {
            NodeID H_id = aux_n[2] + j;
            NodeID G_id = X[2].global_id[j];

            phi[H_id] = potential[j];
            membership[G_id].push_back(H_id);
            global_id[H_id] = G_id;
        }
    }

    PRINT("ADDING EDGES TO H");
    bool has_padding = false;
    std::vector<FullEdge> e;
    for (int i = 0; i < n; i++) {
        if (membership[i].size() != 1)
            has_padding = true;

        for (auto edge: graph.getEdgesOf(i))
            for (auto u: membership[i])
                for (auto v: membership[edge.target])
                    e.emplace_back(u, v, edge.weight);
    }
    PRINT("CREATING H");
    if (has_padding)
        PRINT("HAS PADDING");
    else
        PRINT("NO PADDING");

    auto H = Graph(H_n, e);
    PRINT("DONE CREATING");

    PRINT("    RUNNING LAZY DIJKSTRA: H_n = " << H_n << ", e.size() = " << e.size() << ", n = " << n << ", m = " << m);

    std::optional<Distances> optional_H_potential;
    if (config::pad_use_lazy)
        optional_H_potential = bcf::runLazyDijkstra(H, phi, Orientation::OUT, config::pad_rounds * 2 + 1); // A max number of rounds
    else
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

// pad
