#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"


bool isCompetitorValid(ImmersedSimplicialComplex competitor, NetEmbeddedWithPins link) {
    auto pin1 = link.embedding[link.pins.first];
    auto pin2 = link.embedding[link.pins.second];
    auto maybe_line = competitor.containsSegment(pin1, pin2);
    if (!maybe_line) {
        cerr << "The line connecting the pins is not contained in the competitor." << endl;
        return false;
    }
    auto line = *maybe_line;

    if (!competitor.isCollapsible(line)) {
        cerr << "The line contains the points:" << endl;
        for (auto singleton: line.getFaces(0))
            cerr << *singleton.begin() << " ";
        cerr << endl;
        cerr << "The competitor is not collapsible on the line connecting the pins." << endl;
        return false;
    }

    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;

        auto U = link.embedding[d.from];
        auto V = link.embedding[d.to];
        if (!competitor.containsSegment(U, V)) {
            cerr << "The competitor does not contain the boundary edge " << U << ", " << V << "." << endl;
            return false;
        }
    }
    return true;
}


bool containsErrorMessage(const string& log_file, const string& error_string) {
    ifstream infile(log_file);
    // const string& pattern = "X Error of failed request";

    assert(infile.is_open());

    string line;
    while (getline(infile, line)) {
        if (line.find(error_string) != std::string::npos)
            return true;
    }
    return false;
}

void processNet(int net_id, bool run_surface_evolver, vector<double> params) {
    clog << "Processing net " << net_id;
    if (!params.empty()) {
        clog << " with params ";
        for (double x: params)
            clog << x << " ";
    }
    clog << endl;

    NetEmbeddedWithPins link = embedNetWithPins(net_id, params);

    if (!pop_description.contains(net_id)) {
        printInfoForHumanPopping(link);
        return;
    }

    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone =
        (net_id == 1029 and params[0] == 1)
        ? generateCubeCompetitor(link)
        : generatePoppedCone(link,
                             pop_description.at(net_id).first,
                             pop_description.at(net_id).second);

    const string& embedding_name = generateEmbeddingFilename(net_id, params);
    const string& filename_prefix = "data/surface_evolver/" + embedding_name;

    saveSurfaceEvolverFormat(
            net_id, embedding_name, cone, filename_prefix + "_cone.fe", false);

    if (run_surface_evolver) {
        saveSurfaceEvolverFormat(
            net_id, embedding_name, cone, filename_prefix + "_cone.fe", false);
        saveSurfaceEvolverFormat(
            net_id, embedding_name, popped_cone, filename_prefix + "_popped.fe", true);

        string surface_evolver_command = format(
            "evolver -x -Q {}_popped.fe > {}.log 2>&1",
            filename_prefix, filename_prefix);
        clog << "Executing the surface evolver..." << endl;
        int surface_evolver_exit_code = system(surface_evolver_command.c_str());
        // The exit code may be nonzero because of some strange interaction
        // between the surface evolver and the X window system when the
        // graphics window is closed with close_show.
        bool x_error = containsErrorMessage(filename_prefix + ".log", "X Error of failed request");
        bool surface_evolver_error = containsErrorMessage(filename_prefix + ".log", "ERROR");
        assert(surface_evolver_exit_code == 0 or (x_error and !surface_evolver_error));

        string convert_format_command = format(
            "/usr/bin/python3  scripts/prettify_surface_evolver_dump.py {}_competitor.out {}_competitor.txt > /dev/null",
            filename_prefix, filename_prefix);
        int convert_format_exit_code = system(convert_format_command.c_str());
        assert(convert_format_exit_code == 0);
    }

    ifstream competitor_in(filename_prefix + "_competitor.txt");
    assert(competitor_in.is_open());
    ImmersedSimplicialComplex competitor = parseImmersedSimplicialComplex(competitor_in);

    competitor = removeAllCollapsedSegments(competitor);
    competitor = removeAllCollapsedTriangles(competitor);

    clog << "Link length: " << link.totalLen() << ", " << link.totalInterestingLen() << endl;
    clog << "Areas" << endl;
    clog << "Cone: " << cone.X.area() << ", "
         << "Popped cone: " << popped_cone.X.area() << ", "
         << "Competitor: " << competitor.area() << endl;

    bool valid_competitor = isCompetitorValid(competitor, link);
    assert(valid_competitor);
    assert(cone.X.area() > competitor.area() + 1e-3);
}

int main() {
    processNet(827, true, {0.46, 0.3});
    return 0;
    for (int net_id: {1029}) {
        for (double d = 0.6; d < acos(1.0/3) / 2 - 0.01; d += 0.005) {
            cout << d << endl;
            processNet(net_id, true, {0., d});
        }
    }

    return 0;
    // testSegmentDistanceFromLine();
    // testTriangleDistanceFromLine();
    // testIsCollapsible();
    // return 0;

    // 827 is the only 2-parameter family.
    // 920 is expected to be optimal.
    // 1029 is the cube family, which is the union of two 1-parameter families.
    set<int> special_nets = {827, 920, 1029};

    // 158 is an instance of 1029 (which we understand algebraically)
    // 681 is an instance of 1029 (same instance as 158)
    // 1100 is the cubic instance of 1029.
    set<int> cube_special_cases = {158, 681, 1100};

    // set<int> heppes_with_pins_on_vertices = {250, 325, 736};
    for (int net_id: possible_nets) {
        if (net_id != 28)
            continue;
        if (special_nets.count(net_id))
            continue;
        if (cube_special_cases.count(net_id))
            continue;

        processNet(net_id, true, {0.03});
    }
} // 5.7273923
