#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

struct SphericalEmbedding {
    vector<vector<int>> adj;
    int n() const { return adj.size(); }
};

struct PlanarEmbedding : SphericalEmbedding {
    vector<int> boundary;
};

// a and b are two vectors.
// Each of them contains distinct values.
// This function checks whether b is a rotation of a.
bool IsRotation(const vector<int>& a, vector<int>& b) {
    if (a.size() != b.size()) return false;
    if (a.empty()) return true;
    auto it = find(b.begin(), b.end(), a[0]);
    if (it == b.end()) return false;
    std::rotate(b.begin(), it, b.end());
    
    return a == b;
}

// Whether a = perm(b)
bool IsRotation(const vector<int>& a, vector<int> b, const vector<int>& perm) {
    for (int& x: b) x = perm[x];
    return IsRotation(a, b);
}

// Checks whether the permutation sigma is an isomorphism from adj1
// to adj2.
bool IsIsomorphism(const SphericalEmbedding& G1,
                   const SphericalEmbedding& G2,
                   const vector<int>& perm) {
    assert(G1.n() == G2.n());
    assert(G1.n() == ssize(perm));
    
    int n = G1.n();
    for (int v = 0; v < n; v++) {
        if (!IsRotation(G2.adj[perm[v]], G1.adj[v], perm)) return false;
    }
    return true;
}

// Checks whether the permutation sigma is an isomorphism from adj1
// to adj2.
bool IsIsomorphism(const PlanarEmbedding& G1,
                   const PlanarEmbedding& G2,
                   const vector<int>& perm) {
    if (!IsIsomorphism(static_cast<SphericalEmbedding>(G1),
                       static_cast<SphericalEmbedding>(G2),
                       perm)) return false;
    
    return IsRotation(G2.boundary, G1.boundary, perm);
}

vector<int> GeometricGraphVisit(const SphericalEmbedding& G, int u, int v) {
    vector<int> ord;
    
    vector<bool> visited(G.n(), false);
    function<void(int,int)> dfs = [&](int u, int v) {
        if (visited[u]) return;
        visited[u] = true;
        ord.push_back(u);
        
        const vector<int>& adj = G.adj[u];
        auto it = find(adj.begin(), adj.end(), v);
        assert(it != adj.end());
        int offset = it - adj.begin();
        int k = ssize(adj);
        for (int i = 0; i < k; i++) {
            int a = adj[(offset + i) % k];
            dfs(a, u);
        }
    };
    
    return ord;
}

// The graphs G1 and G2 are assumed to be connected.
//
// Returns {} if there is not an isomorphism that maps u1 to u2 and v1 to v2.
// Otherwise it returns the unique isomorphism that maps u1 to u2 and v1 to v2.
vector<int> AreIsomorphic(const SphericalEmbedding& G1, int u1, int v1,
                          const SphericalEmbedding& G2, int u2, int v2) {
    if (G1.n() != G2.n()) return {};
    
    vector<int> ord1 = GeometricGraphVisit(G1, u1, v1);
    vector<int> ord2 = GeometricGraphVisit(G2, u2, v2);
    
    vector<int> perm(G1.n());
    for (int i = 0; i < G1.n(); i++) perm[ord1[i]] = ord2[i];
    
    if (IsIsomorphism(G1, G2, perm)) return perm;
    return {};
}

// The graphs G1 and G2 are assumed to be connected.
// Returns {} if there is not an isomorphism.
// Otherwise, it returns an isomorphism.
vector<int> AreIsomorphic(const SphericalEmbedding& G1,
                          const SphericalEmbedding& G2) {
    if (G1.n() != G2.n()) return {};
    for (int v = 0; v < G2.n(); v++) {
        for (int u: G2.adj[v]) {
            vector<int> perm = AreIsomorphic(G1, 0, G1.adj[0][0], G2, v, u);
            if (!perm.empty()) return perm;
        }
    }
    return {};
}

// The graphs G1 and G2 are assumed to be connected.
// Returns {} if there is not an isomorphism.
// Otherwise, it returns an isomorphism.
vector<int> AreIsomorphic(const PlanarEmbedding& G1,
                          const PlanarEmbedding& G2) {
    if (G1.n() != G2.n()) return {};
    
    int k = G2.boundary.size();
    for (int i = 0; i < k; i++) {
        int v = G2.boundary[i];
        int u = G2.boundary[(i + 1) % k];
        vector<int> perm = AreIsomorphic(G1, G1.boundary[0], G1.boundary[1], 
                                         G2, v, u);
        if (!perm.empty()) return perm;
    }
    return {};
}


typedef unsigned long long ULL;
const int HASH_NUM = 1000;
ULL _hash[HASH_NUM];
void fill_hash() {
    srand(0);
    for (int i = 0; i < HASH_NUM; i++) _hash[i] = rand();
}


// We pass also boundary because it avoids repetition when
// planar embeddings are considered instead.
vector<vector<ULL>> EdgesMatrix(const SphericalEmbedding& G,
                                vector<int> boundary = {}) {
    map<pair<int,int>, int> edge2id_memo;
    auto edge2id = [&](int u, int v) {
        if (!edge2id_memo.count({u, v})) {
            int id = edge2id_memo.size();
            edge2id_memo[{u, v}] = id;
        }
        return edge2id_memo[{u, v}];
    };
    
    int E = 0;
    for (int v = 0; v < G.n(); v++)
        E += G.adj[v].size();
    assert(E % 2 == 0);
    E /= 2;
    
    vector<vector<ULL>> mat(E, vector<ULL>(E, 0));
    for (int v = 0; v < G.n(); v++) {
        const auto& adj = G.adj[v];
        int k = adj.size();
        for (int i = 0; i < k; i++) {
            int u = adj[i];
            mat[edge2id(u, v)][edge2id(u, v)] = _hash[HASH_NUM-1];
            for (int j = 0; j < k; j++) {
                int w = adj[j];
                int diff = j-i;
                if (diff < 0) diff += k;
                mat[edge2id(u, v)][edge2id(v, w)] = _hash[diff];
            }
        }
    }
    
    int k = ssize(boundary);
    for (int i = 0; i < k; i++) {
        int u = boundary[i];
        int v = boundary[(i + 1) % k];
        int w = boundary[(i + 2) % k];
        mat[edge2id(u, v)][edge2id(v, w)] += _hash[HASH_NUM-2];
    }
    
    return mat;
}

vector<vector<ULL>> EdgesMatrix(const PlanarEmbedding& G) {
    return EdgesMatrix(G, G.boundary);
}

// path containing l vertices (possibly 0) from u to v.
// u, v must be on the boundary
PlanarEmbedding AddFace(PlanarEmbedding G, int u, int v, int l) {
    auto& b = G.boundary;
    auto it = find(b.begin(), b.end(), v);
    assert(it != b.end());
    it = next(it);
    if (it == b.end()) it = b.begin();
    rotate(b.begin(), it, b.end());
    // At this point, v is the last element of b.
    
    // Updating adj[u] and adj[v].
    int aft = b[ssize(b)-2];
    auto it_u = find(b.begin(), b.end(), u);
    int bef = *next(it_u);
    
    int new_bef = (l == 0) ? v : G.n();
    int new_aft = (l == 0) ? u : G.n() + l - 1;
    
    G.adj[u].insert(find(G.adj[u].begin(), G.adj[u].end(), bef), new_bef);
    G.adj[v].insert(next(find(G.adj[v].begin(), G.adj[v].end(), aft)), new_aft);
    
    // Updating the boundary b.
    while (b.back() != u) b.pop_back();
    for (int i = 0; i < l; i++) {
        int x = G.n();
        b.push_back(x);
        
        vector<int> vec = {x-1, x+1};
        if (i == 0) vec[0] = u;
        if (i == l-1) vec[1] = v;  
        G.adj.push_back(vec);
    }
    b.push_back(v);
    
    return G;
}

int main() {
    fill_hash();
    
    cout << "Ciao" << endl;
    
    
    PlanarEmbedding X;
    AreIsomorphic(X, 0, 0, X, 0, 0); 
}
