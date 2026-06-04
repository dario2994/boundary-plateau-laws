#ifndef MARKED_COMBINATORIAL_MAPS
#define MARKED_COMBINATORIAL_MAPS

#include "utils.cpp"
#include "combinatorial_maps.cpp"

// Class representing a marked combinatorial map.
// It is a combinatorial map such that all vertices have degree
// either 2 or 3.
// Moreover, for all vertices of degree 2, we know which side is convex and
// which side is concave.
//
// This class is used to represent candidates to be once-derived spherical
// graphs. In this context the convex angles are 120 degrees and the concave
// ones are 240.

class MarkedCombinatorialMap : public CombinatorialMap {
  private:
    // If i belongs to _concave_angles, then the angle between
    // dart(i) and phi(dart(i)) is concave.
    vector<Dart> _concave_angles;

  public:
    MarkedCombinatorialMap(const CombinatorialMap& A): CombinatorialMap(A) {}
    friend MarkedCombinatorialMap parseMarkedCombinatorialMap(istream& in);
    friend MarkedCombinatorialMap symmetrize(const MarkedCombinatorialMap&);

    pair<Dart, Dart> addEdge(Dart bef, Dart aft) = delete;
    pair<Dart, Dart> addVertex(Dart bef, Dart aft) = delete;

    bool isConcaveAngle(const Dart& d) const {
        return find(_concave_angles.begin(), _concave_angles.end(), d)
               != _concave_angles.end();
    }

    int concaveAngleId(const Dart& d) const {
        auto it = find(_concave_angles.begin(), _concave_angles.end(), d);
        assert(it != _concave_angles.end());
        return it - _concave_angles.begin();
    }

    int phiAngle(const Dart& d) const {
        return isConcaveAngle(d) ? 60 : -60;
    }

    int nextAngle(const Dart& d) const {
        return isConcaveAngle(rev(d)) ? 240 : 120;
    }

    pair<Dart, Dart> edgeSubdivision(Dart d) {
        bool d_concave = isConcaveAngle(d);
        bool d_rev_concave = isConcaveAngle(rev(d));

        _concave_angles.erase(
            remove(_concave_angles.begin(), _concave_angles.end(), d), _concave_angles.end());
        _concave_angles.erase(
            remove(_concave_angles.begin(), _concave_angles.end(), rev(d)), _concave_angles.end());

        pair<Dart,Dart> pp = CombinatorialMap::edgeSubdivision(d);

        // The newly created angle is concave.
        _concave_angles.push_back(pp.first);

        if (d_concave)
            _concave_angles.push_back(pp.second);
        if (d_rev_concave)
            _concave_angles.push_back(rev(pp.first));

        return pp;
    }

    const vector<Dart>& concaveAngles() const {
        return _concave_angles;
    }

    bool isConvexFace(const Dart& d) const {
        vector<Dart> face = CombinatorialMap::face(d);
        for (const Dart& e: face) {
            if (isConcaveAngle(e))
                return false;
        }
        return true;
    }

    int faceArea(const Dart& d) const {
        const vector<Dart>& f = face(d);
        int res = 6 - ssize(f);
        for (const Dart& e: f) {
            if (isConcaveAngle(e))
                res += 2;
        }
        return res;
    }

    // Returns {# sides, # concave angles} of the face of d.
    pair<int,int> faceType(const Dart& d) const {
        const vector<Dart>& f = face(d);
        int concave_num = 0;
        for (const Dart& e: f) concave_num += isConcaveAngle(e);

        return {ssize(f), concave_num};
    }

    // Whether d and e belong to the same face.
    bool sameFace(const Dart& d, const Dart& e) const {
        const vector<Dart>& f = face(d);
        for (const Dart& d1: f)
            if (d1 == e)
                return true;
        return false;
    }


    bool isConvexFace(const vector<Dart>& f) const {
        return isConvexFace(f[0]);
    }

    int faceArea(const vector<Dart>& f) const {
        return faceArea(f[0]);
    }

    pair<int,int> faceType(const vector<Dart>& f) const {
        return faceType(f[0]);
    }

    ULL isomorphismSignature() const {
        const auto& mat = largeMatrixPower(matrixSignature());

        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][next(d).id]^mat[d.id][rev(d).id];
        }

        for (const Dart& d: concaveAngles()) {
            ULL tmp = 0;
            for (Dart e: darts())
                tmp ^= mat[d.id][e.id];
            res += tmp;
        }
        return res;
    }

    void serialize(ostream& out) const {
        CombinatorialMap::serialize(out);

        out << concaveAngles().size() << "\n";
        for (const Dart& d: concaveAngles())
            out << d.id << " ";
        out << "\n";
    }
};

ostream& operator<<(ostream& out, const MarkedCombinatorialMap& A) {
    out << static_cast<const CombinatorialMap&>(A) << "\n";
    out << "Concave angles: " << "\n";

    for (const Dart& d: A.concaveAngles()) {
        out << d << " " << A.phi(d) << "\n";
    }

    return out;
}

// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the MarkedCombinatorialMaps.
// The ids of the darts are used for a concise representation.
vector<int> findIsomorphism(const MarkedCombinatorialMap& A,
                            const MarkedCombinatorialMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};

    if (A.concaveAngles().size() != B.concaveAngles().size())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = findIsomorphism(A, B, A.dart(0), d);
        if (iso.empty())
            continue;

        bool works = true;
        for (Dart e: A.concaveAngles())
            works &= B.isConcaveAngle(B.dart(iso[e.id]));

        if (works)
            return iso;
    }
    return {};
}

MarkedCombinatorialMap parseMarkedCombinatorialMap(istream& in) {
    MarkedCombinatorialMap A(parseCombinatorialMap(in));
    int concave_angles_num;
    in >> concave_angles_num;

    set<int> concave_vertices;
    for (int i = 0; i < concave_angles_num; i++) {
        int id;
        in >> id;
        A._concave_angles.push_back(A.dart(id));

        int v = A.dart(id).to;
        assert(concave_vertices.count(v) == 0);
        concave_vertices.insert(v);
        assert(A.degree(v) == 2);
    }
    return A;
}

// Reads an entire file and parse all the marked combinatorial maps listed.
vector<MarkedCombinatorialMap> parseMarkedCombinatorialMaps(string in_file) {
    ifstream in(in_file);
    assert(in.is_open());
    int cnt;
    in >> cnt;
    vector<MarkedCombinatorialMap> res;
    for (int i = 0; i < cnt; i++) {
        MarkedCombinatorialMap A = parseMarkedCombinatorialMap(in);
        A.checkConsistency();
        res.push_back(A);
    }
    return res;
}

MarkedCombinatorialMap symmetrize(const MarkedCombinatorialMap& A) {
    vector<Dart> sym_concave_angles;
    for (const Dart& d: A._concave_angles)
        sym_concave_angles.push_back(A.rev(A.phi(d)));

    MarkedCombinatorialMap B(symmetrize(static_cast<CombinatorialMap>(A)));
    B._concave_angles = sym_concave_angles;
    return B;
}


#endif
