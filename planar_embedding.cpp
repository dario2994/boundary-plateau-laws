#ifndef boundary_taylor_planar_embedding
#define boundary_taylor_planar_embedding
#include "utils.cpp"
#include "graph.cpp"


// Represenatation of the embedding of a planar graph in the plane.
// The graph G is the graph itself, with edges of each vertex ordered
// clockwise.
// The vector boundary holds, in clockwise order, the exterior vertices in
// the embedding.
struct PlanarEmbedding {
    Graph G;
    vector<int> boundary; // at least two elements

    // Returns a vector which represents the isomorphism class of the embedding.
    // Two embeddings have the same representing vector if and only if one
    // is mapped into the other through a rename of the vertices.
    // The complexity is O(|boundary|*|edges(G)|).
    vector<int> isomorphism_class() {
        vector<int> ans;
        for (int it = 0; it < (int)boundary.size(); it++) {
            vector<int> ord = G.ordered_dfs(boundary[it], get(boundary, it+1));
            minimal_permutation(ord);
            if (ans.empty() or ord < ans) ans = ord;
        }
        return ans;
    }

    // Bound from below on the area of an isogonal spherical embedding whose
    // dual graph is obtained from this one adding external edges and vertices.
    // The area measure is normalized so that the sphere has area 12.
    // - Vertices of degree 2 have area 4.
    // - Vertices of degree 3 have area 3 (in the dual they are triangular faces),
    // - Vertices of degree 4 have area 2,
    // - Vertices of degree 5 have area 1.
    // - Vertices of degree >= 6 cannot be present.
    int area() const {
        vector<int> deg_multiplicity(6, 0);
        for (int v = 0; v < G.n; v++) {
            int deg = G.edges[v].size();
            if (deg > 5) return 1<<30;
            deg_multiplicity[deg]++;
        }
        for (int v: boundary) {
            int deg = G.edges[v].size();
            deg_multiplicity[deg]--;
        }
        int area = 0;
        if (deg_multiplicity[0] + deg_multiplicity[1] > 0) return 1<<30;
        for (int d = 2; d <= 5; d++) area += (6-d) * deg_multiplicity[d];
        return area + (int)boundary.size();
    }
};

// Append an external new triangle to boundary[it], boundary[it+1].
PlanarEmbedding AppendTriangle(PlanarEmbedding PE, int it) {
    int v = PE.G.n;
    int b_size = PE.boundary.size();
    assert(b_size >= 2 and it < b_size);
    
    int u = get(PE.boundary, it);
    int w = get(PE.boundary, it+1);
    PE.G.add_vertex();
    PE.G.add_edge(u, v);
    PE.G.add_edge(v, w);
    PE.boundary.insert(PE.boundary.begin() + it+1, v);
    return PE;
}

// Make an external triangle with vertices boundary[it], boundary[it+1],
// boundary[it+2].
PlanarEmbedding CreateTriangle(PlanarEmbedding PE, int it) {
    int b_size = PE.boundary.size();
    assert(b_size >= 3 and it < b_size);
    
    int u = get(PE.boundary, it);
    int w = get(PE.boundary, it+2);
    PE.G.add_edge(u, w);
    PE.boundary.erase(PE.boundary.begin() + (it+1) % b_size);
    return PE;
}

vector<PlanarEmbedding> generate_triangulations(const int max_faces) {
    vector<vector<PlanarEmbedding>> triangulations(max_faces + 1);

    Graph triangle;
    triangle.add_vertex(), triangle.add_vertex(), triangle.add_vertex();
    triangle.add_edge(0, 1), triangle.add_edge(1, 2), triangle.add_edge(2, 0);
    PlanarEmbedding TriangleEmbedding;
    TriangleEmbedding.G = triangle;
    TriangleEmbedding.boundary = {0, 1, 2};
    triangulations[1].push_back(TriangleEmbedding);

    cout << "Generating triangulations on the plane." << endl;
    cout << "Faces\t\tIsomorphism classes" << endl;
    
    for (int faces = 2; faces <= max_faces; faces++) {
        auto add_candidate = [&](const PlanarEmbedding& PE) {
            if (PE.area() > 12) return;
            triangulations[faces].push_back(PE);
        };
        for (PlanarEmbedding PE: triangulations[faces-1]) {
            int boundary_size = PE.boundary.size();
            for (int i = 0; i < boundary_size; i++) {
                add_candidate(AppendTriangle(PE, i));
                if (boundary_size > 2) add_candidate(CreateTriangle(PE, i));
            }
        }
        remove_isomorphic_copies(triangulations[faces]);

        cout << faces << ":\t\t" << triangulations[faces].size() << endl;
    }
    cout << endl;

    vector<PlanarEmbedding> all_triangulations;
    for (int faces = 1; faces <= max_faces; faces++)
        for (PlanarEmbedding PE: triangulations[faces])
            all_triangulations.push_back(PE);
    
    return all_triangulations;
}

#endif
