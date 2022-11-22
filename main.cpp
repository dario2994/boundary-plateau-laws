#include "utils.cpp"
#include "graph.cpp"
#include "planar_embedding.cpp"
#include "spherical_embedding.cpp"
#include "identify_equal_sides.cpp"

double triangle_side() {
    return acos(-1.0 / 3);
}

// Given the length of one side of a rectangle, returns the length of an
// adjacent side.
// Decreasing function.
// 0 < a < pi. (a = 0 is the area 
// The limiting case a = 0 (b = pi) is the section delimited by two meridians
// forming an angle of 60 degres.
double rectangle_side(double a) {
    double cosb = (3-5*cos(a)) / (5-3*cos(a));
    return acos(cosb);
}

// Given the lengths of two adjacent sides of a pentagon, returns the length
// of side adjacent to neither.
// Increasing function (with respect to both parameters).
double pentagon_side(double a, double b) {
    double cosd = (1.0/3 +  cos(a) + cos(b) + cos(a) * cos(b) - 2*sin(a)*sin(b)) / 2;
    return acos(cosd);
}

// If it returns true, then there is at most one isogonal embedding up to
// isometries.
bool unique_isogonal_embedding(const SphericalEmbedding& SE) {
    map<pair<int,int>, int> label = equal_sides(SE);
    vector<bool> is_known(label.size(), false);

    vector<vector<int>> faces;
    for (int v = 0; v < SE.G.n; v++) {
        vector<int> face;
        for (int u: SE.G.edges[v]) face.push_back(label[{v, u}]);
        faces.push_back(face);
    }

    while (1) {
        bool changed = false;
        bool finished = true;
        for (const vector<int>& face: faces) {
            bool is_regular = true;
            int cnt_notknown = 0;
            for (int x: face) {
                cnt_notknown += !is_known[x];
                is_regular &= x == face[0];
            }
            if (cnt_notknown == 0) continue;
            finished = false;
            if (cnt_notknown <= 3 or is_regular) {
                changed = true;
                for (int x: face) is_known[x] = true;
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
    double z = acos(sqrt(5)/3);
    cout << pentagon_side(z, z) << " " << z << endl;
    cout << triangle_side() << endl;
    
    for (double a = 0.1; a < 3; a += 0.2) {
        double from = 10;
        double to = 0;
        for (double b = 0.01; b < 6; b += 0.001) {
            double d = pentagon_side(a, b);
            if (!isnan(d)) {
                from = min(from, b);
                to = max(to, b);
            }
        }
        for (double b = 0.1; b < 3; b += 0.1) {
            double d = pentagon_side(a, b);
            if (from < b and b < to) assert(!isnan(d));
        }
        cout << a << ": " << from << " " << to << endl;
        cout << a + from << endl;
    }
    // cout << rectangle_side(1) << endl;
    return 0;
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
