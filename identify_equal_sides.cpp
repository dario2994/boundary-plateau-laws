#ifndef boundary_taylor_identify_equal_sides
#define boundary_taylor_identify_equal_sides
#include "spherical_embedding.cpp"

bool find_contradiction_assuming_ineq(
        vector<vector<int>> faces, vector<pair<int,int>> ineqs,
        vector<int> label) {
    map<int,int> label2id;
    auto get_id = [&](int t) {
        if (label2id.count(t) == 0) {
            int sz = label2id.size();
            label2id[t] = sz;
        }
        return label2id[t];
    };
    for (int& t: label) t = get_id(t);
    for (auto& f: faces) for (int& t: f) t = label[t];
    
    int q = label2id.size();
    vector<vector<bool>> is_less(q, vector<bool>(q, false));
    for (auto ineq: ineqs)
        is_less[get_id(ineq.first)][get_id(ineq.second)] = true;

    // From now on we care only about faces and is_less.

    // This function completes the is_less relation by applying the implication
    // is_less[i][j] and is_less[j][k] implies is_less[i][k].
    auto transitive_closure = [&]() {
        bool changed = false;
        for (int i = 0; i < q; i++) for (int j = 0; j < q; j++) {
            if (!is_less[i][j]) continue;
            for (int k = 0; k < q; k++) {
                if (is_less[j][k] and !is_less[i][k]) {
                    is_less[i][k] = true;
                    changed = true;
                }
            }
        }
        return changed;
    };

    // Comparing pentagons.
    auto heppes_ineq_lemma = [&](const vector<int>& a, const vector<int>& b) {
        if (a.size() != 5 or b.size() != 5) return false;
        if (a[0] != b[0]) return false;

        int parity = -1;
        for (int i = 1; i < 5; i++) {
            int p = a[i], q = b[i];
            if (is_less[p][q]) parity = i%2;
        }
        if (parity == -1) return false;

        bool changed = false;
        for (int i = 1; i < 5; i++) {
            int p = a[i], q = b[i];
            if (i%2 != parity) swap(p, q);
            if (!is_less[p][q]) {
                changed = true;
                is_less[p][q] = true;
            }
        }
        return changed;
    };

    while (1) {
        bool changed = false;
        for (const auto& a: faces) for (const auto& b: faces)
            changed |= heppes_ineq_lemma(a, b);
        changed |= transitive_closure();

        if (!changed) break;
    }

    for (int i = 0; i < q; i++) for (int j = 0; j < q; j++)
        if (is_less[i][j] and is_less[j][i]) return true; // contradiction
    
    return false;
}

// Returns a map edge -> label, such that if two edges have the same label
// then they have the same length for any isogonal spherical embedding SE.
map<pair<int,int>, int> equal_sides(const SphericalEmbedding& SE) {
    const auto& G = SE.G;
    
    map<pair<int,int>, int> edge2id;
    int e_num = 0;
    for (int v = 0; v < G.n; v++) {
        for (int u: G.edges[v]) {
            if (u < v) continue;
            edge2id[{v, u}] = edge2id[{u, v}] = e_num++;
        }
    }

    vector<vector<int>> faces;
    for (int v = 0; v < G.n; v++) {
        vector<int> x;
        for (int u: G.edges[v]) x.push_back(edge2id[{v, u}]);
        int q = x.size();
        for (int direction: {1, -1}) for (int pos = 0; pos < q; pos++) {
            vector<int> y(q);
            for (int i = 0; i < q; i++) y[i] = get(x, pos + direction * i);
            faces.push_back(y);
        }
    }

    vector<int> label(e_num);
    for (int i = 0; i < e_num; i++) label[i] = i;
    auto join = [&](int x, int y) {
        x = label[x], y = label[y];
        for (int i = 0; i < e_num; i++) if (label[i] == x) label[i] = y;
    };

    auto prove_isometric = [&](const vector<int>& x, const vector<int>& y) {
        assert(x.size() == y.size());
        int q = x.size();

        int cnt_nonequal = 0;
        for (int i = 0; i < q; i++) cnt_nonequal += label[x[i]] != label[y[i]];

        if (0 < cnt_nonequal and cnt_nonequal <= 3) {
            for (int i = 0; i < q; i++) join(x[i], y[i]);
            return true;
        }
        return false;
    };

    bool is_dodecahedron =
        SE.faces_number() == 12 and SE.edges_number() == (12*5)/2;

    if (is_dodecahedron) {
        // Consider two adjacent faces in the dodecahedron
        //   _ _
        //  <_|_>
        // 
        // which correspond to two adjacent vertices in the dual.
        // Assume that the central edge is strictly larger than all
        // the other ones.
        // If we find a contradiction, then we may assume that there's
        // a symmetric pentagon (which is equivalent to assuming that
        // there are two adjacent edges with the same length).
        vector<pair<int,int>> ineqs;
        int u = 0, v = G.edges[u][0];
        
        int id = edge2id[{u, v}];
        for (int a: G.edges[u]) 
            if (a != v) ineqs.push_back({id, edge2id[{u, a}]});
        for (int a: G.edges[v]) 
            if (a != u) ineqs.push_back({id, edge2id[{v, a}]});
        
        if (find_contradiction_assuming_ineq(faces, ineqs, label))
            join(id, edge2id[{u, G.edges[u][1]}]);
    }
        
    while (1) {
        bool changed = false;
        for (const auto& x: faces) for (const auto& y: faces) {
            if (x.size() != y.size()) continue;
            changed |= prove_isometric(x, y);
        }
        if (changed) continue;

        vector<int> label_vals = label;
        sort(label_vals.begin(), label_vals.end());
        label_vals.erase(unique(label_vals.begin(), label_vals.end()),
                         label_vals.end());

        for (int x: label_vals) for (int y: label_vals) {
            if (y <= x) continue;
            if (find_contradiction_assuming_ineq(faces, {{x, y}}, label)) {
                changed = true;
                join(x, y);
                break;
            }
            if (changed) break;
        }
        
        if (!changed) break;
    }

    for (auto& pp: edge2id) pp.second = label[pp.second];
    return edge2id;
}

#endif
