#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"


ErrorOr<bool> isCompetitorValid(ImmersedSimplicialComplex competitor, NetEmbeddedWithPins link) {
    auto pin1 = link.embedding[link.pins.first];
    auto pin2 = link.embedding[link.pins.second];
    auto maybe_line = competitor.containsSegment(pin1, pin2);
    if (!maybe_line.isOk()) {
        return ErrorMessage("The line connecting the pins is not properly contained in the competitor.\n"
                            + maybe_line.getError());
    }
    auto line = maybe_line.getContent();

    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;

        auto U = link.embedding[d.from];
        auto V = link.embedding[d.to];
        auto maybe_boundary_edge = competitor.containsSegment(U, V);
        if (!maybe_boundary_edge.isOk()) {
            stringstream ss;
            ss << "The competitor does not contain the boundary edge " << U << ", " << V << "." << endl;
            return ErrorMessage(ss.str() + maybe_line.getError());
        }
    }

    auto maybe_collapsible = competitor.isCollapsible(line);
    if (!maybe_collapsible.isOk()) {
        stringstream ss;
        ss << "The line contains the points:" << endl;
        for (auto singleton: line.getFaces(0))
            ss << *singleton.begin() << " ";
        ss << endl;
        ss << "The competitor is not collapsible on the line connecting the pins." << endl;
        return ErrorMessage(ss.str() + maybe_collapsible.getError());
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

SurfaceEvolverData generatePoppedConeCases(int net_id, vector<double> params, NetEmbeddedWithPins link) {
    if (net_id == 1029 and params[0] == 1) // Cube family
        return generateCubeCompetitor(link);
    if (net_id == 827) {// Fat fish family // TODO: Bad!
        double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(params[1])));
        if (params[0] > upperbound_theta0 - 0.05)
            return generateFatFishCompetitor2(link);
        else
            return generateFatFishCompetitor(link);
    }
    return generatePoppedCone(link, pop_description.at(net_id).first,
                                    pop_description.at(net_id).second);
}

ErrorOr<double> processNet(int net_id, bool run_surface_evolver, vector<double> params) {
    NetEmbeddedWithPins link = embedNetWithPins(net_id, params);

    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedConeCases(net_id, params, link);

    const string& embedding_name = generateEmbeddingFilename(net_id, params);
    const string& filename_prefix = "data/surface_evolver/" + embedding_name;

    if (run_surface_evolver) {
        saveSurfaceEvolverFormat(
            net_id, embedding_name, cone, filename_prefix + "_cone.fe", false);
        saveSurfaceEvolverFormat(
            net_id, embedding_name, popped_cone, filename_prefix + "_popped.fe", true);

        string surface_evolver_command = format(
            "evolver -x -Q {}_popped.fe > {}.log 2>&1",
            filename_prefix, filename_prefix);
        // clog << "Executing the surface evolver..." << endl;
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

    // competitor = removeAllCollapsedSegments(competitor); // Maybe still necessary?
    // competitor = removeAllCollapsedTriangles(competitor); // Maybe still necessary?

    // clog << "Link length: " << link.totalLen() << ", " << link.totalInterestingLen() << endl;
    // clog << "Areas" << endl;
    // clog << "Cone: " << cone.X.area() << ", "
    //      << "Popped cone: " << popped_cone.X.area() << ", "
    //      << "Competitor: " << competitor.area() << endl;

    CHECK_OR_RETURN(isCompetitorValid(competitor, link));

    double area_gap = cone.X.area() - competitor.area();
    if (area_gap < 1e-3) {
        stringstream ss;
        ss << "The gap = " << area_gap << " between the area of the cone = "
           << cone.X.area() << " and the area of the competitor = "
           << competitor.area() << " is not sufficiently large.";
        return ErrorMessage(ss);
    }

    return area_gap;
}

bool proveInstabilityOneParameterNet(int net_id) {
    Interval param_region;
    if (rigid_nets_flexible_pins.count(net_id)) {
        ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
        assert(metrized_net_in.is_open());
        MetrizedNet net = parseMetrizedNet(metrized_net_in);
        param_region = validPinEdgeLength(net_id, net);
    } else if (net_id == 1029) // goggles
        param_region = {0., acos(1.0/3) / 2};
    else assert(0);

    clog << "Param region: " << param_region << endl;

    auto get_area_gap = [&](double param) {
        if (net_id == 1029) // TODO: Bad!
            return VALUE_OR_ASSERT(processNet(net_id, true, {0., param}));
        if (net_id == 1096 and param > 0.38) { // TODO: Bad!
            pop_description[1096].first[0] = 8;
            double area_gap = VALUE_OR_ASSERT(processNet(net_id, true, {param}));
            pop_description[1096].first[0] = 6;
            return area_gap;
        }
        return VALUE_OR_ASSERT(processNet(net_id, true, {param}));
    };

    auto width_covered = [&](double param) {
        double area_gap = get_area_gap(param);
        clog << "area_gap(" << param << ") = " << area_gap << endl;
        return area_gap / 10.;
    };

    // We execute run_surface_evolver(param) for some values of param
    // in param_region in order to prove the instability for all the possible
    // values of param.
    // Our methods fail if param is too close to one of the extreme points
    // of the interval of admissible values, so we try to avoid such situation
    // as much as possible (this is the reason why the code is a bit involved).
    double new_l  = param_region.mid();
    while (width_covered(new_l) < new_l - param_region.l)
        new_l = (new_l + param_region.l) / 2;

    double new_r  = param_region.mid();
    while (width_covered(new_r) < param_region.r - new_r)
        new_r = (new_r + param_region.r) / 2;

    clog << "new_l = " << new_l << ", new_r = " << new_r << endl;

    // We have proven instability in [I.l, new_l] and [new_r, I.r].
    // So, now we can restrict our analysis to the interval [new_l, new_r].
    double param = new_l;
    while (param < new_r) {
        param += width_covered(param);
    }
    return true;
}

double uniformSample(double x) {
    static std::mt19937 gen(42);

    std::uniform_real_distribution<double> dist(0.0, x);
    return dist(gen);
}


int main() {
    const double sigma = sqrt(2. / 3);
    // const double tau = 1. / sqrt(3.);

    for (int _ = 0; _ < 100; _++) {
        double upperbound_theta1 = acos(sigma);
        double theta1 = uniformSample(upperbound_theta1);
        double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1)));
        double theta0 = uniformSample(upperbound_theta0);

        cout << _ << endl;
        cout << theta0 << ", " << theta1 << endl;
        cout << upperbound_theta0 << " " << upperbound_theta1 << endl;

        auto maybe_proof = processNet(827, true, {theta0, theta1});
        if (!maybe_proof.isOk()) {
            cerr << maybe_proof.getError() << endl;
            assert(0);
        }
    }
    // testSegmentDistanceFromLine();
    // testTriangleDistanceFromLine();
    // testIsCollapsible();
    // return 0;
}