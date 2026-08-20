#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"


struct FindCompetitorConfig {
    int net827_pop_type; // -1, 0, 1, fat fish
    int net1096_removed_face; // Default is 6, sometimes overwritten as 8
    bool net827_slow_convergence; // whether we want to run for a longer time the surface evolver, fat fish
    // TODO
    // bool net1029_cube_family; // Whether we should treat the net 1029 as the cube family. The net 1029 admits two 1-parameter families of metrization, one of them is the cube family, the other is the goggles.
};


// TODO: We should write down the correct competitor and then recheck without changing it.
ErrorOr<bool> isCompetitorValid(GeodesicNetAntipodalPins link, ImmersedSimplicialComplex competitor) {
    vector<Point> embedding = link.embedding;
    Point pin1 = embedding[link.pins.first];
    Point pin2 = embedding[link.pins.second];
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

        Point U = embedding[d.from];
        Point V = embedding[d.to];
        auto maybe_boundary_edge = competitor.containsSegment(U, V, false);
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
    assert(infile.is_open());

    string line;
    while (getline(infile, line)) {
        if (line.find(error_string) != std::string::npos)
            return true;
    }
    return false;
}

SurfaceEvolverData generatePoppedConeCases(int net_id, const GeodesicNetAntipodalPins& link, FindCompetitorConfig config) {
    // if (net_id == 1029 and config.net1029_cube_family)
    //     return generateCubeCompetitor(link); // TODO
    if (net_id == 1096)
        return generatePoppedCone(link, {config.net1096_removed_face}, pop_description.at(net_id).second);
    if (net_id == 827) {// Fat fish family
        return generateFatFishCompetitor(link, config.net827_pop_type);
    }
    return generatePoppedCone(link, pop_description.at(net_id).first,
                                    pop_description.at(net_id).second);
}


ErrorOr<double> findCompetitor(int net_id, string name,
                               const GeodesicNetAntipodalPins& link,
                               bool run_surface_evolver,
                               FindCompetitorConfig config) {
    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedConeCases(net_id, link, config);

    const string& net_directory = format("data/surface_evolver/{}", net_id);
    const string& filename_prefix = format("{}/{}", net_directory, name);

    if (run_surface_evolver) {
        filesystem::create_directory(net_directory);
        saveSurfaceEvolverFormat(cone, filename_prefix + "_cone.fe", "");
        string surface_evolver_commands = surfaceEvolverCommands(
            net_id, filename_prefix, config.net827_slow_convergence);
        saveSurfaceEvolverFormat(
            popped_cone, filename_prefix + "_popped.fe", surface_evolver_commands);

        string bash_command = format(
            "evolver -x -Q {}_popped.fe > {}.log 2>&1",
            filename_prefix, filename_prefix);
        int surface_evolver_exit_code = system(bash_command.c_str());
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

    auto maybe_validation_error = isCompetitorValid(link, competitor);
    if (!maybe_validation_error.isOk())
        return ErrorMessage(maybe_validation_error.getError());

    // TODO: This shall be fixed!
    double area_gap = cone.X.area() - competitor.area();
    if (area_gap < 1e-5) {
        stringstream ss;
        ss << "The gap = " << area_gap << " between the area of the cone = "
           << cone.X.area() << " and the area of the competitor = "
           << competitor.area() << " is not sufficiently large. ";
        return ErrorMessage(ss);
    }

    return area_gap;
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
ErrorOr<double> excludeFamilyOfNets(int net_id, vector<Interval> net_params,
                                    bool run_surface_evolver,
                                    FindCompetitorConfig config) {
    NetEmbeddedWithPins family = embedNetWithPins(net_id, net_params, false);
    GeodesicNetAntipodalPins link =
            embedNetWithPins(net_id, net_params, true).deterministicEmbedding();
    const string& name = generateEmbeddingFilename(net_id, net_params);

    ErrorOr<double> res = findCompetitor(net_id, name, link, run_surface_evolver, config);

    if (!res.isOk())
        return res;

    double area_gap = res.getContent();
    double needed_gap = excludeFamilyViaStability(link, family);

    cout.precision(7);
    cout << fixed << area_gap << "\t(" << needed_gap << ")" << endl;

    assert(needed_gap > 1e-8);
    if (area_gap < needed_gap) {
        stringstream ss;
        ss << "The area gap = " << area_gap << " found for the selected link "
           << "is not sufficient to exclude the family. "
           << "A gap of " << needed_gap << " would be required.";
        return ErrorMessage(ss);
    }

    return area_gap;
}

// If this function terminates, it is guaranteed that there is a finite set family
// of subintervals I_1, .., I_k such that:
// - They cover I, that is I \subseteq \bigcup_j I_j
// - It holds f(I_j) = true for all j.
void coverInterval(Interval I, function<bool(Interval)> f) {
    if (f(I))
        return;
    assert(I.len() > 1e-5); // TODO: Avoid going too deep.
    Interval left = {I.l, I.mid()};
    Interval right = {I.mid(), I.r};
    left.intersect(I);
    right.intersect(I);
    coverInterval(left, f);
    coverInterval(right, f);
};

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
        auto maybe_proof = excludeFamilyOfNets(827, {Interval(theta0), Interval(theta1)}, true, config);
        if (!maybe_proof.isOk()) {
            error_messages.push_back(maybe_proof.getError());
        } else {
            area_gap = maybe_proof.getContent();
            if (area_gap > 0.01)
                break;
        }
    }
    // if (ssize(error_messages) == 6) {
    //     for (string msg: error_messages)
    //         cout << msg << endl << endl;
    //     assert(0);
    // }
    cout << " " << area_gap << " (" << ssize(error_messages) << ")";

    return area_gap;
}

// theta0 controls the position of the pins
// theta1 controls the derived-once net
double proveInstabilityFatFishRecursive(Interval theta0, Interval theta1) {
    static double area_processed = 0;

    if (theta0.len() + theta1.len() < 1e-4)
        return -1.0;

    double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1.r)));
    theta0.intersect(Interval(-1000, upperbound_theta0 - 0.01));
    if (theta0.empty())
        return true;
    cout << theta0 << " " << theta1 << ", " << (100*area_processed / 0.36) << "%, ";

    double desired_area_gap =
        // 8 * theta0.len() + 2*sqrt(3) * theta1.len() + 1e-3;
        10 * theta0.len() + 4 * theta1.len() + 1e-4;
    double p0 = theta0.mid();
    double p1 = theta1.mid();

    if (theta0.len() + theta1.len() < 0.01) {
        double area_gap = processFatFish(p0, p1);
        if (desired_area_gap / 2 < area_gap) {
            cout << " OK" << endl;
            area_processed += theta0.len() * theta1.len();
            return area_gap;
        }
    }
    cout << endl;

    bool children_are_solved = true;
    double res = -1;
    for (Interval theta0_child: vector<Interval>{{theta0.l, p0}, {p0, theta0.r}}) {
        for (Interval theta1_child: vector<Interval>{{theta1.l, p1}, {p1, theta1.r}}) {
            double child_res = proveInstabilityFatFishRecursive(theta0_child, theta1_child);
            res = max(res, child_res);
            children_are_solved &= child_res > 0;
            if (desired_area_gap < child_res)
                return child_res;
        }
    }
    if (children_are_solved)
        return res;
    return -1.0;
}

bool proveInstabilityFatFish() {
    const double sigma = sqrt(2. / 3);
    double upperbound_theta1 = acos(sigma);
    return proveInstabilityFatFishRecursive(
        Interval(0, 1000),
        Interval(0, upperbound_theta1-0.01)) > 0;
}


void proveInstabilityRigidNets() {
    cout << "Processing rigid nets." << endl;
    cout << "net_id: area gap" << endl;
    for (int net_id: rigid_nets) {
        // if (net_id ! 325)
        //     continue; // TODO: Remove
        cout << net_id << ": \t";
        cout.flush();
        // TODO: Fix
        auto maybe_proof = excludeFamilyOfNets(net_id, {}, false, FindCompetitorConfig());
        if (!maybe_proof.isOk())
            cout << maybe_proof.getError() << endl;
        assert(maybe_proof.isOk());
        // double area_gap = maybe_proof.getContent();
        cout << endl;
    }
    cout << endl;
}

void proveInstabilityRigidNetsWithFlexiblePins() {
    cout << "Processing rigid nets with flexible pins." << endl;
    for (int net_id: rigid_nets_flexible_pins) {
        cout << net_id << endl;

        ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
        assert(metrized_net_in.is_open());
        ApproximatelyMetrizedNet net = parseApproximatelyMetrizedNet(metrized_net_in);
        Interval param_region = validPinEdgeLength(net);
        if (net_id == 1096) {
            param_region.l += 0.005;
            param_region.r -= 0.005; // TODO
        }
        clog << "Param region: " << param_region << endl;

        // Returns true or false.
        auto try_interval = [&](Interval param) {
            if (param.len() > 0.01) // No matter what, we split at least this much.
                return false;
            clog << "Processing " << param << endl;
            FindCompetitorConfig config;
            config.net1096_removed_face = (net_id == 1096 and param.mid() > 0.38) ? 8 : 6;
            return excludeFamilyOfNets(net_id, {param}, true, config).isOk();
        };

        coverInterval(param_region, try_interval);

        cout << "OK" << endl << endl;
    }
}

void proveInstabilityGoggles() {
    cout << "Processing goggles." << endl;
    Interval param_region = {0., acos(1.0/3) / 2 - 0.005}; // TODO: Justify
    clog << "Param region: " << param_region << endl;

    // Returns true or false.
    auto try_interval = [&](Interval param) {
        if (param.len() > 0.01) // No matter what, we split at least this much.
            return false;
        clog << "Processing " << param << endl;
        return excludeFamilyOfNets(GOGGLES_ID, {param}, true, FindCompetitorConfig()).isOk();
    };

    coverInterval(param_region, try_interval);

    cout << "OK" << endl << endl;
}

int main() {
    // proveInstabilityRigidNets();
    // proveInstabilityRigidNetsWithFlexiblePins();
    // proveInstabilityGoggles();
    proveInstabilityFatFish();

    return 0;
}