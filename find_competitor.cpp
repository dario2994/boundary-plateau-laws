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
    if (!maybe_line)
        return false;
    auto line = *maybe_line;

    if (!competitor.isCollapsible(line))
        return false;

    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;

        auto U = link.embedding[d.from];
        auto V = link.embedding[d.to];
        if (!competitor.containsSegment(U, V))
            return false;
    }
    return true;
}


bool containsXServerError(const string& log_file) {
    ifstream infile(log_file);
    const string& pattern = "X Error of failed request";

    assert(infile.is_open());

    string line;
    while (getline(infile, line)) {
        if (line.find(pattern) != std::string::npos)
            return true;
    }
    return false;
}

void processNet(int net_id, bool run_surface_evolver, vector<double> params) {
    NetEmbeddedWithPins link = embedNetWithPins(net_id, params);

    if (!pop_description.contains(net_id)) {
        printInfoForHumanPopping(link);
        return;
    }

    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedCone(link,
        pop_description.at(net_id).first,
        pop_description.at(net_id).second);

    const string& embedding_name = generateEmbeddingFilename(net_id, params);
    const string& filename_prefix = "data/surface_evolver/" + embedding_name;

    if (run_surface_evolver) {
        saveSurfaceEvolverFormat(
            embedding_name, cone, filename_prefix + "_cone.fe", false);
        saveSurfaceEvolverFormat(
            embedding_name, popped_cone, filename_prefix + "_popped.fe", true);

        string surface_evolver_command = format(
            "evolver -x -Q {}_popped.fe > {}.log 2>&1",
            filename_prefix, filename_prefix);
        cout << "Executing the surface evolver..." << endl;
        int surface_evolver_exit_code = system(surface_evolver_command.c_str());
        // The exit code may be nonzero because of some strange interaction
        // between the surface evolver and the X window system when the
        // graphics window is closed with close_show.
        bool x_error = containsXServerError(filename_prefix + ".log");
        assert(surface_evolver_exit_code == 0 or x_error);

        string convert_format_command = format(
            "/usr/bin/python3  scripts/prettify_surface_evolver_dump.py {}_competitor.out {}_competitor.txt > /dev/null",
            filename_prefix, filename_prefix);
        int convert_format_exit_code = system(convert_format_command.c_str());
        assert(convert_format_exit_code == 0);
    }

    ifstream competitor_in(filename_prefix + "_competitor.txt");
    ImmersedSimplicialComplex competitor = parseImmersedSimplicialComplex(competitor_in);

    // cout << "#vertices" << endl;
    // cout << "Cone: " << cone.X.V() << ", "
    //      << "Popped cone: " << popped_cone.X.V() << ", "
    //      << "Competitor: " << competitor.V() << ", ";
    // cout << endl;

    cout << "Areas" << endl;
    cout << "Cone: " << cone.X.area() << ", "
         << "Popped cone: " << popped_cone.X.area() << ", "
         << "Competitor: " << competitor.area() << ", ";

    bool valid_competitor = isCompetitorValid(competitor, link);
    assert(valid_competitor);
    assert(cone.X.area() > competitor.area() + 1e-3);
}

int main() {
    // testSegmentDistanceFromLine();
    // testTriangleDistanceFromLine();
    // testIsCollapsible();
    // return 0;

    processNet(182, true, {});
}
