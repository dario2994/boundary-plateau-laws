#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"

struct FindCompetitorConfig {
    int net827_pop_type; // -1, 0, 1, fat fish
    int net1096_removed_face; // Default is 6, sometimes overwritten as 8
    bool net827_slow_convergence; // whether we want to run for a longer time the surface evolver, fat fish
    bool net1029_cube_family; // Whether we should treat the net 1029 as the cube family. The net 1029 admits two 1-parameter families of metrization, one of them is the cube family, the other is the goggles.
};


// TODO: We should write down the correct competitor and then recheck without changing it.
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
            return ErrorMessage(ss.str() + maybe_boundary_edge.getError());
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

SurfaceEvolverData generatePoppedConeCases(int net_id, NetEmbeddedWithPins link, FindCompetitorConfig config) {
    if (net_id == 1029 and config.net1029_cube_family)
        return generateCubeCompetitor(link); // TODO
    if (net_id == 1096)
        return generatePoppedCone(link, {config.net1096_removed_face}, pop_description.at(net_id).second);
    if (net_id == 827) {// Fat fish family
        return generateFatFishCompetitor(link, config.net827_pop_type);
    }
    return generatePoppedCone(link, pop_description.at(net_id).first,
                                    pop_description.at(net_id).second);
}

// Tries to find a valid competitor for the given net.
// The net is specified via net_id and, if the net is not rigid, via some real
// parameters net_params.
// One can decide to skip running the surface evolver and used the cached
// competitors. Even in this case, this function checks whether the cached
// competitors are valid.
// If the function finds a valid competitor, it saves it and returns the (positive)
// difference between the cone and the competitor.
// Otherwise, it returns an error.
ErrorOr<double> processNet(int net_id, vector<double> net_params,
                           bool run_surface_evolver,
                           FindCompetitorConfig config) {
    NetEmbeddedWithPins link = embedNetWithPins(net_id, net_params);

    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedConeCases(net_id, link, config);

    const string& embedding_name = generateEmbeddingFilename(net_id, net_params);
    const string& filename_prefix = "data/surface_evolver/" + embedding_name;

    if (run_surface_evolver) {
        saveSurfaceEvolverFormat(
            net_id, embedding_name, cone, filename_prefix + "_cone.fe", false, config.net827_slow_convergence);
        saveSurfaceEvolverFormat(
            net_id, embedding_name, popped_cone, filename_prefix + "_popped.fe", true, config.net827_slow_convergence);

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
        if (net_id == 1029) {
            FindCompetitorConfig config;
            config.net1029_cube_family = false;
            return VALUE_OR_ASSERT(processNet(net_id, {0., param}, true, config));
        }
        if (net_id == 1096 and param > 0.38) {
            FindCompetitorConfig config;
            config.net1096_removed_face = (param > 0.38) ? 8 : 6;
            return VALUE_OR_ASSERT(processNet(net_id, {param}, true, config));
        }
        return VALUE_OR_ASSERT(processNet(net_id, {param}, true, FindCompetitorConfig()));
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

double processFatFish(double theta0, double theta1) {
    const static vector<tuple<int, bool>> config_options = {
        {0, false}, {1, false}, {1, true}, {0, true}, {-1, false}, {-1, true}
    };

    double area_gap = -1;
    vector<string> error_messages;
    for (auto config_data: config_options) {
        FindCompetitorConfig config;
        config.net827_pop_type = get<0>(config_data);
        config.net827_slow_convergence = get<1>(config_data);
        auto maybe_proof = processNet(827, {theta0, theta1}, true, config);
        if (!maybe_proof.isOk()) {
            error_messages.push_back(maybe_proof.getError());
        } else {
            area_gap = maybe_proof.getContent();
            if (area_gap > 0.01)
                break;
        }
    }
    if (ssize(error_messages) == 6) {
        for (string msg: error_messages)
            cout << msg << endl << endl;
        assert(0);
    }
    cout << " " << area_gap << " (" << ssize(error_messages) << ")";

    return area_gap;
}

double area_processed = 0;
bool proveInstabilityFatFishRecursive(Interval theta0, Interval theta1) {
    double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1.r)));
    theta0.intersect(Interval(-1000, upperbound_theta0 - 0.01));
    if (theta0.empty())
        return true;
    cout << theta0 << " " << theta1 << ", " << (100*area_processed / 0.36) << "%, ";

    double p0 = theta0.mid();
    double p1 = theta1.mid();

    if (theta0.len() + theta1.len() < 0.01) {
        double area_gap = processFatFish(p0, p1);
        if ((theta0.len() + theta1.len()) / 2 < area_gap / 10) {
            cout << endl;
            area_processed += theta0.len() * theta1.len();
            return true;
        }
    }
    cout << endl;

    auto res00 = proveInstabilityFatFishRecursive({theta0.l, p0}, {theta1.l, p1});
    auto res01 = proveInstabilityFatFishRecursive({theta0.l, p0}, {p1, theta1.r});
    auto res10 = proveInstabilityFatFishRecursive({p0, theta0.r}, {theta1.l, p1});
    auto res11 = proveInstabilityFatFishRecursive({p0, theta0.r}, {p1, theta1.r});

    return res00 and res01 and res10 and res11;
}

bool proveInstabilityFatFish() {
    const double sigma = sqrt(2. / 3);
    double upperbound_theta1 = acos(sigma);
    // double upperbound_theta0 = atan(1.0 / (sqrt(3.)));
    return proveInstabilityFatFishRecursive(
        Interval(0, 1000),
        Interval(0, upperbound_theta1-0.01));
}


double uniformSample(double x) {
    static std::mt19937 gen(42);

    std::uniform_real_distribution<double> dist(0.0, x);
    return dist(gen);
}

void proveRigidNets() {
    set<int> nets_to_ignore = {
        158, // The same net as 1029.
        681, // The same net as 1029.
        920, // Thin fish (proven algebraically)
        1100, // A cube (proven algebraically).
    };

    cout << "net_id: area gap" << endl;
    for (int net_id: possible_nets) {
        if (net_id == 827) // fat fish
            continue;
        if (net_id == 1029) // goggles
            continue;
        if (rigid_nets_flexible_pins.count(net_id))
            continue;

        if (nets_to_ignore.count(net_id))
            continue;

        cout << net_id << ": ";
        cout.flush();
        auto maybe_proof = processNet(net_id, {}, true, FindCompetitorConfig());
        assert(maybe_proof.isOk());
        double area_gap = maybe_proof.getContent();
        cout << area_gap << endl;
    }
}

void proveNetWithOneFreeParameter() {
    vector<int> nets_to_process(rigid_nets_flexible_pins.begin(),
                                rigid_nets_flexible_pins.end());
    nets_to_process.push_back(1029);
    for (int net_id: nets_to_process) {
        cout << net_id << ": ";
        cout.flush();
        assert(proveInstabilityOneParameterNet(net_id));
        cout << "OK" << endl << endl;
    }
}


int main() {
    proveRigidNets();

    proveNetWithOneFreeParameter();

    proveInstabilityFatFish();
    return 0;
}