#ifndef boundary_taylor_spherical_embedding
#define boundary_taylor_spherical_embedding
#include "utils.cpp"
#include "graph.cpp"
#include "planar_embedding.cpp"

struct SphericalEmbedding {
    Graph G; // Dual graph!

    // See SphericalEmbedding::isomorphism_class() for a description.
    vector<int> isomorphism_class() {
        vector<int> ans;
        for (int u = 0; u < G.n; u++) for (int v: G.edges[u]) {
            vector<int> ord = G.ordered_dfs(u, v);
            minimal_permutation(ord);
            if (ans.empty() or ord < ans) ans = ord;
        }
        return ans;
    }

    // Total area of any isogonal embedding (of the dual of G).
    // The area measure is normalized so that the whole sphere has area 12.
    // See PlanarEmbedding::Area for the formula.
    int area() {
        int ans = 0;
        for (int v = 0; v < G.n; v++) {
            int deg = G.degree(v);
            if (deg < 2 or deg > 5) return 1<<30;
            ans += 6-deg;
        }
        return ans;
    }

    // These utility function are not for G, but for its dual (that is what
    // we really want to embed in the sphere).
    int faces_number() const {
        return G.n;
    }

    int edges_number() const {
        int ans = 0;
        for (int v = 0; v < G.n; v++) ans += G.degree(v);
        return ans / 2;
    }

    int vertices_number() const {
        return edges_number() + 2 - faces_number(); // Euler's formula
    }


    // Exclusion Criteria.
    // If any of this criteria is not met, then the graph does not admit
    // an isogonal embedding.
    
    // Heppes proves that a triangle and a pentagon cannot be both present (in
    // the dual graph).
    bool if_triangle_then_no_pentagon() {
        for (int u = 0; u < G.n; u++) for (int v = 0; v < G.n; v++) {
            if (G.degree(u) == 3 and G.degree(v) == 5) return false;
        }
        return true;
    }

    // If there is a 2-face, then there can be only 2-faces (in the dual graph).
    bool if_twoface_then_fixed() {
        for (int u = 0; u < G.n; u++) for (int v = 0; v < G.n; v++) {
            if (G.degree(u) == 2 and G.degree(v) != 2) return false;
        }
        return true;
    }
};

vector<SphericalEmbedding>
generate_spherical_triangulations(const int max_faces) {
    vector<PlanarEmbedding> triangulations =
            generate_triangulations(max_faces-1);
    
    vector<SphericalEmbedding> spherical_triangulations;

    cout << "Generating triangulations on the sphere." << endl;

    for (PlanarEmbedding PE: triangulations) {
        if (PE.boundary.size() != 3) continue;
        SphericalEmbedding SE;
        SE.G = PE.G;
        if (SE.area() == 12
                and SE.if_triangle_then_no_pentagon()
                and SE.if_twoface_then_fixed())
            spherical_triangulations.push_back(SE);
    }
    remove_isomorphic_copies(spherical_triangulations);

    cout << "Found " << spherical_triangulations.size()
         << " isomorphism classes." << endl << endl;
    
    return spherical_triangulations;
}

#endif
