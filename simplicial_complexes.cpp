#ifndef SIMPLICIAL_COMPLEXES
#define SIMPLICIAL_COMPLEXES

#include "utils.cpp"
#include "interval_geometry.cpp"

void assertBoundedNorm(const vector<Point>& points) {
    for (const Point& P: points)
        assert(norm(P) < 3);
}

void assertWellDefinedLine(const Point& A, const Point& B) {
    assertBoundedNorm({A, B});
    assert(norm(A-B) > 1e-4);
}

// Distance of the point P to the line through A and B.
// Assumption: |A-B| is not much smaller than max(|A|, |B|).
double distanceFromLine(Point A, Point B, Point P) {
    assertWellDefinedLine(A, B);
    assertBoundedNorm({P});
    Point L = (B-A) / norm(B-A);
    P = P - A;

    return norm(L ^ P);
}

// Distance of the segment CD to the line through A and B.
// Assumption: |A-B| is not much smaller than max(|A|, |B|).
double distanceFromLine(const Point& A,
                        const Point& B,
                        const Point& C,
                        const Point& D) {
    assertWellDefinedLine(A, B);
    assertBoundedNorm({C, D});

    Point x = (B - A) ^ (D - C);
    Point y = (B - A) ^ (A - C);

    double s_unclamped = (norm(x) < 1e-5) ? 0.0 : (x * y) / (x * x);
    double s = max(0.0, min(1.0, s_unclamped));

    Point closest_point = C + (D - C) * s;

    return distanceFromLine(A, B, closest_point);
}

// Distance of the triangle CDE to the line through A and B.
// Assumption: |A-B| is not much smaller than max(|A|, |B|).
double distanceFromLine(const Point& A,
                        const Point& B,
                        const Point& C,
                        const Point& D,
                        const Point& E) {
    assertWellDefinedLine(A, B);
    assertBoundedNorm({C, D, E});

    double boundary_distance = min({
        distanceFromLine(A, B, C, D),
        distanceFromLine(A, B, D, E),
        distanceFromLine(A, B, E, C),
    });

    // The answer is either boundary_distance or it is 0.
    // Determining whether is zero is annoying.
    // We use the Moller-Trumbore algorithm.
    // The implementation is inspired from wikipedia.
    Point CD = D - C;
    Point CE = E - C;
    Point CA = A - C;

    Point U = (B-A) ^ CE;
    Point V = CA ^ CD;

    double det = CD * U;

    // There is a plane that "almost" contains all points.
    if (abs(det) < 1e-5)
        return boundary_distance;

    double inv_det = 1.0 / det;
    double u = inv_det * (CA * U);
    double v = inv_det * ((B-A) * V);

    // The intersection between the line and the plane containing the triangle
    // falls outside of the triangle.
    if (u <= 0 or v <= 0 or 1 <= u + v)
        return boundary_distance;

    return 0.0;
}


// Class representing an abstract simplicial complex with at least one vertex.
// In essence, this is a family of subsets of {0, 1, ..., V-1} which contains
// all singletons and is closed with respect to subset.
// On top of this, every vertex is mapped to a point in R^3. By linearity, this
// induces an embedding of the topological simplicial complex in R^3.
struct ImmersedSimplicialComplex {
  private:
    // The vertices are indexed 0, 1, ..., V-1.
    vector<Point> points;
    // All faces, not only maximal ones.
    // Any constructor will make sure that faces are distinct and form a family
    // closed with respect to subset.
    vector<set<int>> faces;

    void generateSubsetClosure() {
        set<set<int>> new_faces;
        for (set<int> face_set: faces) {
            int k = ssize(face_set);
            vector<int> face_vec(face_set.begin(), face_set.end());
            for (int bb = 0; bb < (1<<k); bb++) {
                set<int> subface;
                for (int i = 0; i < k; i++) {
                    if (bb&(1<<i))
                        subface.insert(face_vec[i]);
                }
                new_faces.insert(subface);
            }
        }
        faces = vector<set<int>>(new_faces.begin(), new_faces.end());
    }

  public:
    int V() const { return ssize(points); }

    const Point& getPoint(int id) const { return points[id]; }

    // Get the skeleton of the given dimension.
    // If no dimension is specified, all faces are returned.
    vector<set<int>> getFaces(int dim = -1234) const {
        vector<set<int>> res;
        for (const set<int>& face: faces) {
            if (dim == -1234 or ssize(face) == dim+1)
                res.push_back(face);
        }
        return res;
    }

    ImmersedSimplicialComplex(int _V, vector<set<int>> faces): points(_V), faces(faces) {
        for (const auto& face: faces)
            for (int v: face)
                assert(0 <= v and v < V());
        generateSubsetClosure();
    }
    ImmersedSimplicialComplex(vector<array<double,3>> _points,
                              vector<set<int>> faces): ImmersedSimplicialComplex(ssize(_points), faces) {
        points = _points;
    }

    friend ostream& operator<<(ostream&, const ImmersedSimplicialComplex&);

    double area() const {
        double total_area = 0;
        for (auto face: faces) {
            if (ssize(face) != 3)
                continue;
            vector<Point> vertices;
            for (int v: face)
                vertices.push_back(points[v]);
            array<double,3> triangle_vec = (vertices[2]-vertices[0])^(vertices[1]-vertices[0]);
            double triangle_area = 0.5 * norm(triangle_vec);
            total_area += triangle_area;
        }
        return total_area;
    }

    // Whether the simplex can be collapsed on the target subcomplex (on the same vertices).
    ErrorOr<bool> isCollapsible(const ImmersedSimplicialComplex& target) const {
        int F = ssize(faces);
        map<set<int>, int> face2id;
        for (const auto& face: faces)
            face2id[face] = ssize(face2id);

        vector<bool> is_target(F, false);
        for (const set<int>& target_face: target.faces)
            is_target[face2id[target_face]] = true;

        vector<vector<int>> boundary(F);
        vector<vector<int>> coboundary(F);
        for (const auto& face: faces) {
            int id = face2id[face];
            for (int v: face) {
                set<int> subface = face;
                subface.erase(v);
                assert(face2id.count(subface));
                int sub_id = face2id[subface];
                boundary[id].push_back(sub_id);
                coboundary[sub_id].push_back(id);
            }
        }
        queue<int> free_faces;
        vector<bool> removed(F, false);
        for (int id = 0; id < F; id++) {
            if (ssize(coboundary[id]) == 1 and !is_target[id])
                free_faces.push(id);
        }

        auto remove_face = [&](int id) {
            assert(!removed[id]);
            assert(coboundary[id].empty());
            removed[id] = true;
            for (int sub_id: boundary[id]) {
                vector<int>& vec = coboundary[sub_id];
                vec.erase(remove(vec.begin(), vec.end(), id), vec.end());
                if (ssize(vec) == 1 and !is_target[sub_id])
                    free_faces.push(sub_id);
            }
        };

        while (!free_faces.empty()) {
            int id = free_faces.front();
            free_faces.pop();
            if (removed[id] or coboundary[id].empty())
                continue;
            assert(ssize(coboundary[id]) == 1);
            int super_id = coboundary[id][0];
            remove_face(super_id);
            remove_face(id);
        }

        vector<set<int>> non_collapsible_faces;

        for (int i = 0; i < F; i++) {
            if (is_target[i])
                assert(!removed[i]);
            else if (!removed[i])
                non_collapsible_faces.push_back(faces[i]);
        }

        set<int> non_collapsible_vertices;
        for (auto face: non_collapsible_faces)
            for (int v: face)
                non_collapsible_vertices.insert(v);

        if (!non_collapsible_faces.empty()) {
            stringstream ss;
            ss << "Non collapsed vertices: ";
            for (int v: non_collapsible_vertices)
                ss << v << " ";
            ss << endl;
            ss << "Non collapsed faces:" << endl;
            for (auto face: non_collapsible_faces) {
                for (int v: face)
                    ss << v << " ";
                ss << endl;
            }
            return ErrorMessage(ss);
        }

        return true;
    }

    ErrorOr<bool> isCollapsible() const {
        return isCollapsible(ImmersedSimplicialComplex(V(), {}));
    }

    // Whether there is a sequence of vertices v_1, v_2, ..., v_k
    // in the simplicial complex such that (we identify the vertices with their
    // embedding):
    // 1. A = v_1
    // 2. B = v_k
    // 3. The sequence belongs, in this order, to the segment AB
    // 4. For all 1 <= i < k, the segment v_iv_{i+1} is a face of the complex.
    // 5. Any face of the complex is either quantitatively "far" from the
    //    segment AB or is a sub/superset of {v_i, v_{i+1}} for some i.
    // 6. There is no triangle with all three vertices belonging to
    //    {v_1, v_2, ..., v_k}.
    //
    // If allow_changes = true, then this function is allowed to modify
    // the complex to satisfy the above assumptions (but sometimes it may fail
    // in any case).
    ErrorOr<ImmersedSimplicialComplex> containsSegment(
                         const Point& A,
                         const Point& B,
                         bool allow_changes=true,
                         const double eps=1e-6) {
        Point segment_dir = (B-A) * (1.0 / norm(B-A));

        // The vertices are mapped on the line AB.
        vector<bool> on_the_line(V(), false);
        vector<int> points_on_line;
        for (int i = 0; i < V(); i++) {
            double dist = distanceFromLine(A, B, points[i]);
            if (dist > 100*eps)
                continue;

            // We move the points so that they are "perfectly" on the line.
            if (allow_changes) {
                points[i] = A + segment_dir * (segment_dir * (points[i] - A));
                dist = distanceFromLine(A, B, points[i]);
            }

            if (dist > eps) {
                stringstream ss;
                ss << "The relationship between the point " << i << " (that is "
                   << points[i] << ") and the line " << A << ", " << B << " is ambiguous. "
                   << "The distance from the point to the line is " << dist << " and the "
                   << "point would be considered on the line if the distance were < "
                   << eps << ".";
                return ErrorMessage(ss);
            }
            points_on_line.push_back(i);
            on_the_line[i] = true;
        }

        // Sorting the vertices by the order in which they appear on AB.
        auto position_on_line = [&](int i) {
            return segment_dir * (points[i] - A);
        };

        sort(points_on_line.begin(), points_on_line.end(), [&](int i, int j) {
            return position_on_line(i) < position_on_line(j);
        });

        if (ssize(points_on_line) < 2)
            return ErrorMessage("There are less than 2 points on the line.");

        // The first vertex is mapped onto A.
        if (norm(points[points_on_line[0]]-A) > eps) {
            stringstream ss;
            ss << "The intersection between the complex and the line does not start at point"
               << " A = " << A << ".";
            return ErrorMessage(ss);
        }

        // The last vertex is mapped onto B.
        if (norm(points[points_on_line.back()]-B) > eps) {
            stringstream ss;
            ss << "The intersection between the complex and the line does not end at point"
                 << " B = " << B << ".";
            return ErrorMessage(ss);
        }

        // We check that A and B are univocally matched with vertices for the
        // simplicial complex.
        for (int i = 0; i < ssize(points_on_line); i++) {
            Point P = points[points_on_line[i]];
            if (i != 0 and norm(P-A) < 100*eps) {
                stringstream ss;
                ss << "There are multiple vertices of the complex which "
                     << "may coincide with A = " << A << ".";
                return ErrorMessage(ss);
            }
            if (i != ssize(points_on_line)-1 and norm(P-B) < 100*eps) {
                stringstream ss;
                ss << "There are multiple vertices of the complex which "
                     << "may coincide with B = " << B << ".";
                return ErrorMessage(ss);
            }
        }

        // At this point, we have identified the vertices on the segment
        // and their order.
        // We need to "clean" the simplicial complex so that the line is
        // nicely contained. This may, a priori, change the topology.
        // Hopefully, in the relevant cases, it does not.
        if (allow_changes) {
            vector<int> vertex2ord(V(), -1);
            for (int i = 0; i < ssize(points_on_line); i++)
                vertex2ord[points_on_line[i]] = i;
            auto is_on_line = [&](int v) { return vertex2ord[v] >= 0; };
            vector<set<int>> new_faces;

            // Segments connecting consecutive vertices on the line.
            for (int i = 0; i < ssize(points_on_line)-1; i++)
                new_faces.push_back({points_on_line[i], points_on_line[i+1]});

            for (set<int> face: faces) {
                if (ssize(face) <= 1)
                    new_faces.push_back(face);
                if (ssize(face) == 2) {
                    int a = *face.begin(), b = *face.rbegin();
                    if (!is_on_line(a) or !is_on_line(b))
                        new_faces.push_back(face);
                }
                if (ssize(face) == 3) {
                    vector<int> face_vec(face.begin(), face.end());
                    vector<int> iol(3);
                    for (int i = 0; i < 3; i++)
                        iol[i] = is_on_line(face_vec[i]);
                    if (iol[0] + iol[1] + iol[2] <= 1)
                        new_faces.push_back(face);
                    if (iol[0] + iol[1] + iol[2] == 2) {
                        int v = -1, a = -1, b = -1;
                        if (iol[0] == 0)
                            v = 0, a = 1, b = 2;
                        if (iol[1] == 0)
                            v = 1, a = 0, b = 2;
                        if (iol[2] == 0)
                            v = 2, a = 0, b = 1;

                        v = face_vec[v];
                        a = vertex2ord[face_vec[a]];
                        b = vertex2ord[face_vec[b]];
                        if (a > b)
                            swap(a, b);

                        for (int i = a; i < b; i++)
                            new_faces.push_back({v, points_on_line[i], points_on_line[i+1]});
                    }
                }
            }
            faces = new_faces;
            generateSubsetClosure();
        }

        // We save the edges of the simplicial complex in a convenient
        // datastructure to speed up the computation.
        set<pair<int,int>> edges;
        for (auto face: faces) {
            if (ssize(face) != 2)
                continue;
            int x = *face.begin();
            int y = *face.rbegin();
            edges.insert({x, y});
            edges.insert({y, x});
        }

        for (int i = 0; i < ssize(points_on_line)-1; i++) {
            int p = points_on_line[i];
            int q = points_on_line[i+1];
            if (!edges.count({p, q})) {
                stringstream ss;
                ss   << "points[" << p << "] = " << points[p] << ", "
                     << "points[" << q << "] = " << points[q]
                     << " are consecutive on the (spatial) segment "
                     << "but are not connected by a (simplicial) segment.";
                return ErrorMessage(ss);
            }
        }

        // Any face that has only vertices in the sequence is either a
        // singleton or an edge of the path.

        stringstream error_message;
        bool works = true;
        for (int i = 0; i < ssize(points_on_line); i++) {
            for (int j = i+2; j < ssize(points_on_line); j++) {
                int u = points_on_line[i];
                int v = points_on_line[j];
                if (edges.count({u, v})) {
                    error_message   << "points[" << u << "] = " << points[u] << ", "
                         << "points[" << v << "] = " << points[v]
                         << " are not consecutive on the (spatial) segment "
                         << "but they are connected by a (simplicial) segment. "
                         << "Their positions on the (spatial) segment are: "
                         << i << " " << j << "." << endl;
                    error_message << "The points in between are: ";
                    for (int k = i+1; k < j; k++) {
                        int w = points_on_line[k];
                        error_message << "points[" << w << "] = " << points[w] << ", ";
                    }
                    error_message << endl;
                    works = false;
                }
            }
        }
        if (!works)
            return ErrorMessage(error_message);

        // All faces (edges and triangles) that don't have any vertex in the
        // sequence do not intersect the line.
        for (auto face: faces) {
            if (ssize(face) <= 1)
                continue;
            bool has_vertex_on_line = false;
            for (int v: face) {
                if (on_the_line[v])
                    has_vertex_on_line = true;
            }
            if (has_vertex_on_line)
                continue;

            assert(ssize(face) <= 3);
            if (ssize(face) == 2) {
                int c = *face.begin(), d = *face.rbegin();
                Point C = points[c];
                Point D = points[d];
                if (distanceFromLine(A, B, C, D) < 100*eps) {
                    stringstream ss;
                    ss   << "points[" << c << "] = " << C << ", "
                         << "points[" << d << "] = " << D
                         << " are not on the line, but they form a segment "
                         << "that intersects the line." << endl;
                    return ErrorMessage(ss);
                }
            }
            if (ssize(face) == 3) {
                vector<int> ids(face.begin(), face.end());
                if (distanceFromLine(A, B,
                                     points[ids[0]], points[ids[1]], points[ids[2]]) < 100*eps) {
                    stringstream ss;
                    ss   << "points[" << ids[0] << "] = " << points[ids[0]] << ", "
                         << "points[" << ids[1] << "] = " << points[ids[1]] << ", "
                         << "points[" << ids[2] << "] = " << points[ids[2]] << ", "
                         << " are not on the line, but they form a triangle "
                         << "that intersects the line." << endl;
                    return ErrorMessage(ss);
                }
            }
        }

        vector<set<int>> segment_edges;
        for (int i = 0; i < ssize(points_on_line)-1; i++)
            segment_edges.push_back({points_on_line[i], points_on_line[i+1]});


        return ImmersedSimplicialComplex(points, segment_edges);
    }
};

ImmersedSimplicialComplex parseImmersedSimplicialComplex(istream& in) {
    int V, F;
    in >> V >> F;
    vector<Point> points(V);
    for (int i = 0; i < V; i++)
        in >> points[i][0] >> points[i][1] >> points[i][2];
    vector<set<int>> faces;
    for (int i = 0; i < F; i++) {
        int k;
        in >> k;
        set<int> face;
        for (int j = 0; j < k; j++) {
            int v;
            in >> v;
            face.insert(v);
        }
        faces.push_back(face);
    }

    return ImmersedSimplicialComplex(points, faces);
}


ostream& operator<<(ostream& out, const ImmersedSimplicialComplex& X) {
    out << "V = " << X.V() << ", faces = " << (ssize(X.faces)-X.V()-1) << endl;
    for (auto face: X.faces) {
        if (ssize(face) <= 1)
            continue;
        for (int v: face) cout << v << " ";
        out << endl;
    }
    return out;
}


#endif
