#include "utils.cpp"
#include "graph.cpp"
#include "planar_embedding.cpp"
#include "spherical_embedding.cpp"
#include "identify_equal_sides.cpp"


// If it returns true, then there is at most one isogonal embedding up to
// isometries.
// TODO: Can be cleaned using faces as in the previous function.
bool unique_isogonal_embedding(const SphericalEmbedding& SE) {
    map<pair<int,int>, int> label = equal_sides(SE);
    vector<bool> is_known(label.size(), false);
    auto make_known = [&](int u, int v) {
        is_known[label[{u, v}]] = true;
        is_known[label[{v, u}]] = true;
    };

    while (1) {
        bool changed = false;
        bool finished = true;
        for (int v = 0; v < SE.G.n; v++) {
            bool is_regular = true;
            int cnt_notknown = 0;
            for (int u: SE.G.edges[v]) {
                int l = label[{v, u}];
                cnt_notknown += !is_known[l];
                is_regular &= l == label[{v, SE.G.edges[v][0]}];
            }
            if (cnt_notknown == 0) continue;
            finished = false;
            if (cnt_notknown <= 3 or is_regular) {
                changed = true;
                for (int u: SE.G.edges[v]) make_known(u, v);
            }
        }
        if (finished) return true;
        if (!changed) break;
        // TODO: Here we might compute the degrees of freedom by making an
        // arbitrary label known.
    }
    return false;
}

int main() {
    const int max_faces = 20;
    vector<SphericalEmbedding> spherical_triangulations =
            generate_spherical_triangulations(max_faces);
 
    int cnt = 0;
    for (SphericalEmbedding SE: spherical_triangulations) {
        cout << cnt++ << "\t";
        cout << "Vertices: " << SE.vertices_number() << ",\t";
        cout << "Edges: " << SE.edges_number() << ",\t";
        cout << "Faces: " << SE.faces_number() << ",\t";
        cout << "Unique isogonal embedding: "
             << unique_isogonal_embedding(SE) << ",\t";
        set<int> sides;
        auto labels = equal_sides(SE);
        for (auto pp: labels) sides.insert(pp.second);
        cout << "Different lengths: " << sides.size() << ",\t";
        cout << endl;
        
        cout << "Degrees: ";
        for (int v = 0; v < SE.G.n; v++) cout << SE.G.degree(v) << " ";
        cout << endl;

        // SE.G.pretty_print();
    }
}

// TODO: Prove the rigidity of the two embeddings which have two choices of the
//       side lengths due to symmetries. This is easy with some algebra, but
//       may be messy to implement (if one wants also a proof).
//       Create pictures of the 9 cases.
//       Improve the comments, in particular in identify_equal_sides.cpp
//       Make this a git repository urgently.
//
//
// Properties of isogonal embeddings that we are using:
// 1. Apart from the triplane, there are only triangles, quadrilateral, pentagons.
// 2. Triangles are metrically determined.
// 3. Rectangles are metrically determined as soon as one side is fixed.
// 4. If there's a triangle, then cannot be a pentagon (follows from Heppes' Lemma 1).
// 5. The area of a faces with s sides is (6-s)/12 the area of the whole sphere.
// 
// Heppes' Lemmata:
// 1. The side of any pentagon is strictly smaller than the side of any triangle.
// 2. Let (a_i)_{i=1,2,3,4,5} and (b_i)_{i=1,2,3,4,5} be the ordered sides of
//    two pentagons.
//    If a_1 = b_1, we have three cases:
//      - a_i = b_i for all i, so the pentagons are equal.
//      - a_2 < b_2, a_3 > b_3, a_4 < b_4, a_5 > b_5
//      - a_2 > b_2, a_3 < b_3, a_4 > b_4, a_5 < b_5
//
// 
