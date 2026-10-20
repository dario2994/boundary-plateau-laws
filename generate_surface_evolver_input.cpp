#ifndef GENERATE_SURFACE_EVOLVER_INPUT
#define GENERATE_SURFACE_EVOLVER_INPUT


#include "utils.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "simplicial_complexes.cpp"
#include "geodesic_net_antipodal_pins.cpp"
#include "popped_cones.cpp"


string generateEmbeddingFilename(int net_id, vector<Interval> params) {
    std::ostringstream oss;
    oss << net_id;
    oss << fixed << setprecision(5);
    for (Interval param: params)
        oss << "_" << param.l << "_" << param.r;
    return oss.str();
}


string surfaceEvolverCommands(int net_id, string filename_prefix, bool net827_slow_convergence) {
    const string generic_surface_evolver_commands = R"(
    r 2;
    subroutine 5;
    o;
    subroutine 3;
    o;
    subroutine 3;
    r;
    u 4;
    subroutine 5;
    o;
    subroutine 3;)";

    const map<int, string> special_surface_evolver_commands = {
        {28, R"(
            o;
            subroutine 5;
        )"},
        {99, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {250, R"(
            r;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
        )"},
        {346, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {421, R"(
            subroutine 3;
            r;
            subroutine 5;
            r;
            subroutine 3;
            o;
            subroutine 8;
            o;
            subroutine 3;
            o;
            subroutine 3;
        )"},
        {520, R"(
            r 2;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
        )"},
        {614, R"(
            subroutine 3;
            r;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            r;
            subroutine 3;
        )"},
        {736, R"(
            r 2;
            subroutine 12;
        )"},
        {779, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {827, R"(
            r;
            g 400;
            u 2;
            g 300;
        )"},
        {1029, R"(
            subroutine 3;
            r;
            subroutine 2;
        )"},
        {1096, R"(
            r; subroutine 5; o; subroutine 5; o; subroutine 5; o; subroutine 5;
        )"},
    };

    string surface_evolver_commands = generic_surface_evolver_commands;
    if (special_surface_evolver_commands.count(net_id))
        surface_evolver_commands = special_surface_evolver_commands.at(net_id);

    if (net_id == 827 and net827_slow_convergence) // fat fish
        surface_evolver_commands += "\n u 2;\n g 100;\n r;\n g 1000;";

    constexpr string_view format_str = R"(read
subroutine := {{
    quiet;
    g 2000;
    areaweed 0.0005;
    u 4;
    print total_area;
    quiet off;
}}

subroutinequick := {{
    quiet;
    g 200;
    areaweed 0.0005;
    u 4;
    print total_area;
    quiet off;
}}

gogo := {{
{}
    dump "{}_competitor.out";
    /*close_show;*/ /* This raises an error if the window is not visible. */
}}

/* showq; */
gogo;

quit;
)";
    return format(format_str, surface_evolver_commands, filename_prefix);
}


void saveSurfaceEvolverFormat(const SurfaceEvolverData& data,
                              string out_file,
                              string surface_evolver_commands) {
    const auto& X = data.X;
    Point line = X.getPoint(data.pins.first);

    map<pair<int,int>, int> edge2id;
    auto get_edge = [&](int a, int b) {
        assert(a != b);
        assert(0 <= a and a < X.V());
        assert(0 <= b and b < X.V());

        if (!edge2id.count({a, b})) {
            int id = 1 + edge2id.size() / 2;
            edge2id[{a, b}] = id;
            edge2id[{b, a}] = -id;
        }
        return edge2id[{a, b}];
    };

    auto is_vertex_on_line = [&](int v) {
        return norm(X.getPoint(v) ^ line) < 0.001;
    };

    vector<array<int, 3>> triangles;
    for (auto triangle: X.getFaces(2)) {
        vector<int> triangle_vec(triangle.begin(), triangle.end());
        auto A = triangle_vec[0], B = triangle_vec[1], C = triangle_vec[2];
        triangles.push_back({get_edge(A, B), get_edge(B, C), get_edge(C, A)});
    }
    vector<bool> is_vertex_fixed(X.V(), false);
    for (auto edge: data.boundary_edges) {
        is_vertex_fixed[edge.first] = true;
        is_vertex_fixed[edge.second] = true;
    }

    ofstream out(out_file);

    out << "boundary 1 parameters 1" << endl;
    out << "x1:  " << line[0] << " * p1" << endl;
    out << "x2:  " << line[1] << " * p1" << endl;
    out << "x3:  " << line[2] << " * p1" << endl;
    out << endl;

    out << "vertices" << endl;
    for (int i = 0; i < X.V(); i++) {
        Point P = X.getPoint(i);
        bool fixed = is_vertex_fixed[i];
        bool boundary = is_vertex_on_line(i);
        boundary &= !fixed;
        out << i+1;
        if (boundary) {
            out << " " << (P * line) << " boundary 1";
        } else {
            out << " " << P[0] << " " << P[1] << " " << P[2];
            if (fixed)
                out << " fixed";
        }
        out << endl;
    }
    out << endl;

    auto is_edge_fixed = [&](pair<int,int> edge) {
        int u = edge.first, v = edge.second;
        return is_vertex_fixed[u] and is_vertex_fixed[v];
    };

    auto is_edge_boundary = [&](pair<int,int> edge) {
        int u = edge.first, v = edge.second;
        return is_vertex_on_line(u) and is_vertex_on_line(v);
    };

    out << "edges" << endl;
    for (auto pp: edge2id) {
        if (pp.second < 0)
            continue;
        out << pp.second << " "
            << pp.first.first+1 << " " << pp.first.second+1;
        if (is_edge_fixed(pp.first))
            out << " fixed";
        if (is_edge_boundary(pp.first))
            out << " boundary 1";
        out << endl;
    }
    out << endl;

    out << "faces" << endl;
    for (int i = 0; i < ssize(triangles); i++) {
        array<int, 3> T = triangles[i];
        out << i+1 << " "
            << T[0] << " " << T[1] << " " << T[2] << endl;
    }
    out << endl;

    out << R"(read
print total_area;
)";

    out << surface_evolver_commands;

    out.close();
}

#endif
