#include <bits/stdc++.h>
using namespace std;

/*
  Deterministic negative-weight SSSP, engineering implementation inspired by
  Li (2026), "Deterministic Padded Decompositions and Negative-Weight Shortest Paths".

  What is implemented:
    - iterative potential scaling: weights >= -B -> compute potential for shifted graph w+B/2
    - recursive padded-decomposition potential routine
    - deterministic ball-growing decomposition with in/out balls and heavy case
    - layered Hybrid Bellman-Ford/Dijkstra for paths with <= eta negative edges
    - final Johnson/Dijkstra pass

  Practical note:
    The paper's asymptotic bound uses Thorup integer priority queues and the exact BNWN
    scaling framework. This standalone implementation uses std::priority_queue and a final
    Bellman-Ford validation/cleanup to make the public solve() robust. It is meant as a
    faithful, efficient research prototype, not a certified near-linear complexity proof.
*/

struct NegativeSSSP {
    using Real = long double;
    static constexpr Real INF = 1e100L;
    static constexpr Real EPS = 1e-12L;

    struct Edge {
        int u, v;
        Real w;
    };

    struct Result {
        bool has_negative_cycle = false;
        vector<Real> dist;
        vector<int> parent;        // predecessor in a shortest-path tree, -1 if none
    };

    int n;
    vector<Edge> edges;
    vector<vector<pair<int,int>>> outIds, inIds; // (to, edge id), (from, edge id)

    explicit NegativeSSSP(int n = 0) : n(n), outIds(n), inIds(n) {}

    void add_edge(int u, int v, long long w) {
        int id = (int)edges.size();
        edges.push_back({u, v, (Real)w});
        outIds[u].push_back({v, id});
        inIds[v].push_back({u, id});
    }

private:
    struct WorkEdge {
        int u, v;
        Real w; // nonnegative decomposition weight or arbitrary H weight
    };

    struct DijkstraOut {
        vector<Real> dist;
        vector<int> parentVertex;
    };

    struct Decomp {
        vector<vector<int>> sets;
        vector<Real> volumes;
    };

    struct BallInfo {
        Real r = 0;
        vector<int> verts;
        Real vol = 0;
    };

    static int ceil_log2_int(int x) {
        int r = 0;
        int v = 1;
        while (v < max(1, x)) { v <<= 1; ++r; }
        return r;
    }

    static Real safe_log(Real x) {
        return max<Real>(1.0L, log(max<Real>(x, 3.0L)));
    }

    // Multi-source Dijkstra over explicitly supplied nonnegative adjacency.
    static DijkstraOut dijkstra_multi(
        int N,
        const vector<vector<pair<int, Real>>>& adj,
        const vector<Real>& initial,
        bool keep_parent = false
    ) {
        vector<Real> dist = initial;
        vector<int> par(N, -1);
        using QN = pair<Real,int>;
        priority_queue<QN, vector<QN>, greater<QN>> pq;
        for (int i = 0; i < N; ++i)
            if (dist[i] < INF/4) pq.push({dist[i], i});
        while (!pq.empty()) {
            auto [du, u] = pq.top(); pq.pop();
            if (du != dist[u]) continue;
            for (auto [v, w] : adj[u]) {
                Real nd = du + max<Real>(0, w);
                if (nd + EPS < dist[v]) {
                    dist[v] = nd;
                    if (keep_parent) par[v] = u;
                    pq.push({nd, v});
                }
            }
        }
        return {move(dist), move(par)};
    }

    static vector<Real> dijkstra_source(
        int N,
        const vector<vector<pair<int, Real>>>& adj,
        int s
    ) {
        vector<Real> init(N, INF);
        init[s] = 0;
        return dijkstra_multi(N, adj, init).dist;
    }

    // Bellman-Ford from an implicit super-source with zero edges to every active vertex.
    // Returns h such that w(u,v)+h[u]-h[v] >= 0 on the active induced subgraph.
    bool bellman_ford_potential(
        const vector<int>& verts,
        const vector<char>& active,
        const vector<Real>& curW,
        Real shift,
        vector<Real>& phi
    ) const {
        phi.assign(n, 0);
        vector<int> idOf(n, -1);
        for (int i = 0; i < (int)verts.size(); ++i) idOf[verts[i]] = i;
        vector<Edge> sub;
        sub.reserve(edges.size());
        for (int ei = 0; ei < (int)edges.size(); ++ei) {
            const auto& e = edges[ei];
            if (active[e.u] && active[e.v]) sub.push_back({e.u, e.v, curW[ei] + shift});
        }
        bool any = false;
        int k = (int)verts.size();
        for (int it = 0; it < k; ++it) {
            any = false;
            for (const auto& e : sub) {
                if (phi[e.u] + e.w + EPS < phi[e.v]) {
                    phi[e.v] = phi[e.u] + e.w;
                    any = true;
                }
            }
            if (!any) break;
        }
        if (any) {
            // One more pass: still relaxable => negative cycle in shifted graph.
            for (const auto& e : sub)
                if (phi[e.u] + e.w + EPS < phi[e.v]) return false;
        }
        return true;
    }

    // Dijkstra ball on an induced subgraph. direction=+1: out-ball, -1: in-ball.
    BallInfo ball_grow(
        const vector<int>& verts,
        const vector<char>& allowed,
        const vector<Real>& vol,
        const vector<WorkEdge>& wEdges,
        const vector<vector<pair<int,int>>>& wOut,
        const vector<vector<pair<int,int>>>& wIn,
        int s,
        int direction,
        Real delta0,
        Real delta,
        Real eps,
        Real M
    ) const {
        Real lg = safe_log(M + 2);
        Real step = max<Real>(delta * eps / lg, EPS);
        Real limit = delta0 + delta + step;

        vector<Real> dist(n, INF);
        priority_queue<pair<Real,int>, vector<pair<Real,int>>, greater<pair<Real,int>>> pq;
        if (!allowed[s]) return {delta0, {}, 0};
        dist[s] = 0;
        pq.push({0, s});
        while (!pq.empty()) {
            auto [du, u] = pq.top(); pq.pop();
            if (du != dist[u]) continue;
            if (du > limit + EPS) continue;
            const auto& adj = (direction > 0 ? wOut[u] : wIn[u]);
            for (auto [to, id] : adj) {
                const WorkEdge& e = wEdges[id];
                int v = to;
                if (!allowed[v]) continue;
                Real nd = du + e.w;
                if (nd + EPS < dist[v] && nd <= limit + EPS) {
                    dist[v] = nd;
                    pq.push({nd, v});
                }
            }
        }

        vector<pair<Real,int>> byDist;
        byDist.reserve(verts.size());
        for (int v : verts) if (allowed[v] && dist[v] < INF/4) byDist.push_back({dist[v], v});
        sort(byDist.begin(), byDist.end());

        auto volume_at = [&](Real r) {
            Real ans = 0;
            for (auto [dv, v] : byDist) {
                if (dv <= r + EPS) ans += vol[v]; else break;
            }
            return ans;
        };
        auto vertices_at = [&](Real r) {
            vector<int> ans;
            for (auto [dv, v] : byDist) {
                if (dv <= r + EPS) ans.push_back(v); else break;
            }
            return ans;
        };

        const Real C = 8.0L;
        int maxI = (int)ceil(lg / max<Real>(eps, 1e-6L)) + 4;
        Real chosen = delta0;
        for (int i = 1; i <= maxI; ++i) {
            Real rPrev = delta0 + (i-1) * step;
            Real rCur  = delta0 + i * step;
            if (rPrev > delta0 + delta + EPS) break;
            Real vp = max<Real>(1.0L, volume_at(rPrev));
            Real vc = volume_at(rCur);
            if (vc <= (1.0L + C * eps) * vp + EPS) {
                chosen = min<Real>(rPrev, delta0 + delta);
                return {chosen, vertices_at(chosen), volume_at(chosen)};
            }
        }
        chosen = delta0 + delta;
        return {chosen, vertices_at(chosen), volume_at(chosen)};
    }

    vector<int> multi_source_ball(
        const vector<int>& verts,
        const vector<char>& active,
        const vector<int>& sources,
        const vector<WorkEdge>& wEdges,
        const vector<vector<pair<int,int>>>& wOut,
        const vector<vector<pair<int,int>>>& wIn,
        int direction,
        Real radius
    ) const {
        vector<Real> dist(n, INF);
        priority_queue<pair<Real,int>, vector<pair<Real,int>>, greater<pair<Real,int>>> pq;
        for (int s : sources) if (active[s]) {
            dist[s] = 0;
            pq.push({0, s});
        }
        while (!pq.empty()) {
            auto [du, u] = pq.top(); pq.pop();
            if (du != dist[u]) continue;
            if (du > radius + EPS) continue;
            const auto& adj = (direction > 0 ? wOut[u] : wIn[u]);
            for (auto [to, id] : adj) {
                int v = to;
                if (!active[v]) continue;
                Real nd = du + wEdges[id].w;
                if (nd + EPS < dist[v] && nd <= radius + EPS) {
                    dist[v] = nd;
                    pq.push({nd, v});
                }
            }
        }
        vector<int> ans;
        for (int v : verts) if (active[v] && dist[v] <= radius + EPS) ans.push_back(v);
        return ans;
    }

    Real set_volume(const vector<int>& set, const vector<Real>& vol) const {
        Real s = 0;
        for (int v : set) s += vol[v];
        return s;
    }

    Decomp padded_decomposition(
        const vector<int>& verts,
        const vector<char>& active,
        const vector<Real>& curW,
        Real shift,
        Real d,
        Real eps
    ) const {
        // Build the nonnegative graph G_{>=0}^{shift}[X].
        vector<WorkEdge> wEdges;
        wEdges.reserve(edges.size());
        vector<vector<pair<int,int>>> wOut(n), wIn(n);
        vector<Real> vol(n, 0);
        for (int ei = 0; ei < (int)edges.size(); ++ei) {
            const auto& e = edges[ei];
            if (!active[e.u] || !active[e.v]) continue;
            Real ww = max<Real>(0, curW[ei] + shift);
            int id = (int)wEdges.size();
            wEdges.push_back({e.u, e.v, ww});
            wOut[e.u].push_back({e.v, id});
            wIn[e.v].push_back({e.u, id});
            vol[e.u] += 1;
            if (e.v != e.u) vol[e.v] += 1;
        }
        for (int v : verts) vol[v] = max<Real>(1.0L, vol[v]); // regularization for isolates
        Real totalVol = 0;
        for (int v : verts) totalVol += vol[v];
        Real M = max<Real>(1.0L, totalVol / 2.0L);
        Real Delta = max<Real>(d / 12.0L, EPS);
        Real pad = max<Real>(eps * Delta / safe_log(M + 2), EPS);

        auto heavy_case = [&](int s) -> Decomp {
            BallInfo bp = ball_grow(verts, active, vol, wEdges, wOut, wIn, s, +1, Delta, Delta, eps, M);
            BallInfo bm = ball_grow(verts, active, vol, wEdges, wOut, wIn, s, -1, Delta, Delta, eps, M);
            vector<int> BpPad = multi_source_ball(verts, active, vector<int>{s}, wEdges, wOut, wIn, +1, bp.r + pad);
            vector<int> Bp    = multi_source_ball(verts, active, vector<int>{s}, wEdges, wOut, wIn, +1, bp.r);
            vector<int> BmPad = multi_source_ball(verts, active, vector<int>{s}, wEdges, wOut, wIn, -1, bm.r + pad);
            vector<int> Bm    = multi_source_ball(verts, active, vector<int>{s}, wEdges, wOut, wIn, -1, bm.r);
            vector<char> inBpPad(n, false), inBp(n, false), inBmPad(n, false), inBm(n, false);
            for (int v : BpPad) inBpPad[v] = true;
            for (int v : Bp) inBp[v] = true;
            for (int v : BmPad) inBmPad[v] = true;
            for (int v : Bm) inBm[v] = true;
            vector<int> X, Y, Z;
            for (int v : verts) {
                if (inBpPad[v] && inBmPad[v]) X.push_back(v);
                if (inBpPad[v] && !inBm[v]) Y.push_back(v);
                if (!inBp[v]) Z.push_back(v);
            }
            vector<vector<int>> sets;
            for (auto& S : {X, Y, Z}) if (!S.empty()) sets.push_back(S);
            vector<Real> vols;
            for (auto& S : sets) vols.push_back(set_volume(S, vol));
            return {sets, vols};
        };

        vector<char> Uplus(n, false), Uminus(n, false);
        Real volPlus = 0, volMinus = 0;
        bool toHeavy = false;
        int heavySource = verts.empty() ? 0 : verts[0];

        int guard = 0;
        while (volPlus + EPS < M / 2.0L && volMinus + EPS < M / 2.0L && guard++ <= (int)verts.size() + 5) {
            int s = -1;
            for (int v : verts) if (active[v] && !Uplus[v] && !Uminus[v]) { s = v; break; }
            if (s == -1) break;

            vector<char> allowPlus = active, allowMinus = active;
            for (int v : verts) {
                if (Uplus[v]) allowPlus[v] = false;
                if (Uminus[v]) allowMinus[v] = false;
            }
            BallInfo bp = ball_grow(verts, allowPlus, vol, wEdges, wOut, wIn, s, +1, 0, Delta, eps, M);
            BallInfo bm = ball_grow(verts, allowMinus, vol, wEdges, wOut, wIn, s, -1, 0, Delta, eps, M);
            bool choosePlus = (bp.vol <= bm.vol);
            BallInfo chosen = choosePlus ? bp : bm;
            if (chosen.vol <= M + EPS) {
                if (choosePlus) {
                    for (int v : chosen.verts) if (!Uplus[v]) { Uplus[v] = true; volPlus += vol[v]; }
                } else {
                    for (int v : chosen.verts) if (!Uminus[v]) { Uminus[v] = true; volMinus += vol[v]; }
                }
            } else {
                toHeavy = true;
                heavySource = s;
                break;
            }
        }
        if (toHeavy) return heavy_case(heavySource);

        bool usePlus = (volPlus >= volMinus);
        vector<int> U;
        U.reserve(verts.size());
        for (int v : verts) if (usePlus ? Uplus[v] : Uminus[v]) U.push_back(v);
        if (U.empty()) {
            // Degenerate case: force a balanced split to guarantee progress.
            vector<int> A, B;
            for (int i = 0; i < (int)verts.size(); ++i) (i < (int)verts.size()/2 ? A : B).push_back(verts[i]);
            vector<vector<int>> sets;
            if (!A.empty()) sets.push_back(A);
            if (!B.empty()) sets.push_back(B);
            vector<Real> vols;
            for (auto& S : sets) vols.push_back(set_volume(S, vol));
            return {sets, vols};
        }
        vector<int> X = multi_source_ball(verts, active, U, wEdges, wOut, wIn, usePlus ? +1 : -1, pad);
        vector<char> inU(n, false);
        for (int v : U) inU[v] = true;
        vector<int> Y;
        for (int v : verts) if (!inU[v]) Y.push_back(v);
        vector<vector<int>> sets;
        if (!X.empty()) sets.push_back(X);
        if (!Y.empty()) sets.push_back(Y);
        vector<Real> vols;
        for (auto& S : sets) vols.push_back(set_volume(S, vol));
        return {sets, vols};
    }

    struct HybridResult {
        bool ok;
        vector<Real> dist;
    };

    static HybridResult hybrid_bf_dijkstra(
        int N,
        int s,
        const vector<WorkEdge>& hEdges,
        int eta
    ) {
        vector<vector<pair<int, Real>>> nonneg(N);
        vector<WorkEdge> neg;
        for (const auto& e : hEdges) {
            if (e.w >= -EPS) nonneg[e.u].push_back({e.v, max<Real>(0, e.w)});
            else neg.push_back(e);
        }
        vector<Real> dist(N, INF);
        dist[s] = 0;
        dist = dijkstra_multi(N, nonneg, dist).dist;
        for (int k = 1; k <= eta; ++k) {
            vector<Real> next = dist;
            for (const auto& e : neg) {
                if (dist[e.u] < INF/4 && dist[e.u] + e.w + EPS < next[e.v])
                    next[e.v] = dist[e.u] + e.w;
            }
            next = dijkstra_multi(N, nonneg, next).dist;
            dist.swap(next);
        }
        for (const auto& e : hEdges) {
            if (dist[e.u] < INF/4 && dist[e.u] + e.w + EPS < dist[e.v])
                return {false, move(dist)};
        }
        return {true, move(dist)};
    }

    bool recursive_potential(
        const vector<int>& verts,
        Real d,
        const vector<Real>& curW,
        Real shift,
        vector<Real>& phi,
        int depth
    ) const {

        cout << "Recursive Potential: \n";
        cout << "d: " << d << '\n';
        cout << "shift: " << shift << '\n';
        cout << "depth: " << depth << '\n';
        cout << "--------------------------------------------------------------------\n";
        cout.flush();

        phi.assign(n, 0);
        if (verts.empty()) return true;
        vector<char> active(n, false);
        for (int v : verts) active[v] = true;
        int mSub = 0;
        for (int ei = 0; ei < (int)edges.size(); ++ei)
            if (active[edges[ei].u] && active[edges[ei].v]) ++mSub;

        int maxDepth = 2 * ceil_log2_int(max(2, n)) + 30;
        const int BF_EDGE_CUTOFF = 256;
        if (mSub <= BF_EDGE_CUTOFF || (int)verts.size() <= 64 || depth > maxDepth || d < shift + EPS) {
            return bellman_ford_potential(verts, active, curW, shift, phi);
        }

        Real eps = 1.0L / safe_log((Real)max(3, (int)edges.size()));
        Decomp dec = padded_decomposition(verts, active, curW, shift, d, eps);
        if (dec.sets.size() <= 1) return bellman_ford_potential(verts, active, curW, shift, phi);

        // Avoid non-progressive decompositions.
        bool meaningful = false;
        for (const auto& S : dec.sets) if ((int)S.size() < (int)verts.size()) meaningful = true;
        if (!meaningful) return bellman_ford_potential(verts, active, curW, shift, phi);

        vector<vector<Real>> childPhi(dec.sets.size());
        vector<vector<char>> inSet(dec.sets.size(), vector<char>(n, false));
        for (int i = 0; i < (int)dec.sets.size(); ++i) {
            for (int v : dec.sets[i]) inSet[i][v] = true;
            Real di = (dec.volumes[i] <= 1.6L * max<Real>(1.0L, (Real)mSub)) ? d : d / 2.0L;
            if (!recursive_potential(dec.sets[i], di, curW, shift, childPhi[i], depth + 1))
                return false;
        }

        // Build auxiliary graph H.
        vector<vector<int>> copyId(dec.sets.size(), vector<int>(n, -1));
        int Hs = 0, HN = 1;
        for (int i = 0; i < (int)dec.sets.size(); ++i)
            for (int v : dec.sets[i]) copyId[i][v] = HN++;

        vector<WorkEdge> hEdges;
        size_t approx = 0;
        for (auto& S : dec.sets) approx += S.size();
        hEdges.reserve(approx + (size_t)mSub * dec.sets.size() * dec.sets.size());

        for (int i = 0; i < (int)dec.sets.size(); ++i) {
            for (int v : dec.sets[i]) {
                hEdges.push_back({Hs, copyId[i][v], -childPhi[i][v]});
            }
        }
        for (int ei = 0; ei < (int)edges.size(); ++ei) {
            const auto& e = edges[ei];
            if (!active[e.u] || !active[e.v]) continue;
            Real base = curW[ei] + shift;
            for (int i = 0; i < (int)dec.sets.size(); ++i) if (inSet[i][e.u]) {
                int cu = copyId[i][e.u];
                for (int j = 0; j < (int)dec.sets.size(); ++j) if (inSet[j][e.v]) {
                    int cv = copyId[j][e.v];
                    Real ww = base + childPhi[i][e.u] - childPhi[j][e.v];
                    hEdges.push_back({cu, cv, ww});
                }
            }
        }

        int lg = ceil_log2_int(max(3, mSub));
        int eta = min<int>(max<int>(4, 100 * lg * lg + 1), 20000); // cap for practicality
        HybridResult hr = hybrid_bf_dijkstra(HN, Hs, hEdges, eta);
        if (!hr.ok) {
            // The paper extracts a negative cycle here. In this implementation we propagate
            // detection and let the top-level report it.
            return false;
        }

        for (int v : verts) {
            Real best = INF;
            for (int i = 0; i < (int)dec.sets.size(); ++i) if (copyId[i][v] != -1) {
                best = min(best, hr.dist[copyId[i][v]] + childPhi[i][v]);
            }
            if (best >= INF/4) best = 0;
            phi[v] = best;
        }
        return true;
    }

    bool final_bellman_ford_cleanup(vector<Real>& globalPot, vector<Real>& curW) const {
        vector<int> all(n);
        iota(all.begin(), all.end(), 0);
        vector<char> active(n, true);
        vector<Real> h;
        if (!bellman_ford_potential(all, active, curW, 0, h)) return false;
        for (int v = 0; v < n; ++v) globalPot[v] += h[v];
        for (int ei = 0; ei < (int)edges.size(); ++ei) {
            const auto& e = edges[ei];
            curW[ei] += h[e.u] - h[e.v];
            if (curW[ei] < 0 && curW[ei] > -1e-9L) curW[ei] = 0;
        }
        return true;
    }

public:
    Result solve(int source) const {
        Result res;
        res.dist.assign(n, INF);
        res.parent.assign(n, -1);
        if (source < 0 || source >= n) return res;

        vector<Real> curW(edges.size());
        Real minW = INF;
        for (int i = 0; i < (int)edges.size(); ++i) {
            curW[i] = edges[i].w;
            minW = min(minW, curW[i]);
        }
        vector<Real> globalPot(n, 0);
        Real B = max<Real>(0, -minW);

        // Iterative scaling. Each phase computes a potential for current weights shifted by B/2.
        int phases = 0;
        int maxPhases = max(1, 4 + (int)ceil(log2((double)max<Real>(B, 1.0L))) + 4 * ceil_log2_int(n + 2));
        while (B > 1e-7L && phases < maxPhases) {
            Real shift = B / 2.0L;
            vector<int> all(n);
            iota(all.begin(), all.end(), 0);
            Real m = max<Real>(1.0L, (Real)edges.size());
            Real d = m * m * shift;
            vector<Real> psi;
            if (!recursive_potential(all, d, curW, shift, psi, 0)) {
                res.has_negative_cycle = true;
                return res;
            }
            for (int v = 0; v < n; ++v) globalPot[v] += psi[v];
            for (int ei = 0; ei < (int)edges.size(); ++ei) {
                const auto& e = edges[ei];
                curW[ei] += psi[e.u] - psi[e.v];
            }
            B = shift;
            ++phases;
        }

        // Robust final cleanup: makes all reduced costs nonnegative if no negative cycle exists.
        if (!final_bellman_ford_cleanup(globalPot, curW)) {
            res.has_negative_cycle = true;
            return res;
        }

        vector<vector<pair<int, Real>>> adj(n);
        for (int ei = 0; ei < (int)edges.size(); ++ei) {
            const auto& e = edges[ei];
            Real ww = curW[ei];
            if (ww < -1e-8L) {
                res.has_negative_cycle = true;
                return res;
            }
            adj[e.u].push_back({e.v, max<Real>(0, ww)});
        }
        vector<Real> init(n, INF);
        init[source] = 0;
        DijkstraOut dj = dijkstra_multi(n, adj, init, true);
        for (int v = 0; v < n; ++v) {
            if (dj.dist[v] >= INF/4) res.dist[v] = INF;
            else res.dist[v] = dj.dist[v] - globalPot[source] + globalPot[v];
            res.parent[v] = dj.parentVertex[v];
        }
        return res;
    }
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 2 || argc > 3) {
        cerr << "Usage: " << argv[0] << " <graph_file> [source]\n";
        cerr << "Graph format: first line n, then edges u v w until EOF.\n";
        return 1;
    }

    const string graph_file = argv[1];
    const int source = (argc == 3) ? stoi(argv[2]) : 0;

    ifstream in(graph_file);
    if (!in.is_open()) {
        cerr << "Could not open graph file: " << graph_file << '\n';
        return 1;
    }

    int n;
    if (!(in >> n)) {
        cerr << "Could not read node count from graph file: " << graph_file << '\n';
        return 1;
    }

    NegativeSSSP solver(n);
    int u, v;
    long long w;
    while (in >> u >> v >> w) {
        solver.add_edge(u, v, w);
    }

    auto ans = solver.solve(source);
    if (ans.has_negative_cycle) {
        cout << "NEGATIVE CYCLE\n";
        return 0;
    }
    cout.setf(ios::fixed); cout << setprecision(0);
    for (int v = 0; v < n; ++v) {
        if (ans.dist[v] > NegativeSSSP::INF / 8) cout << "INF\n";
        else {
            long double x = ans.dist[v];
            long long rounded = llround((double)x);
            if (fabsl(x - (long double)rounded) < 1e-6L) cout << rounded << "\n";
            else cout << setprecision(12) << (double)x << "\n" << setprecision(0);
        }
    }
    return 0;
}