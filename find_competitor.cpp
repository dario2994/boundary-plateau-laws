#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"


struct FindCompetitorConfig {
    int net827_pop_type; // -1, 0, 1, fat fish
    int net1096_removed_face; // Default is 6, sometimes overwritten as 8
    bool net827_slow_convergence; // whether we want to run for a longer time the surface evolver, fat fish
    bool net1029_cube_family; // Whether we should treat the net 1029 as the cube family.
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
    if (net_id == GOGGLES_ID and config.net1029_cube_family)
        return generateCubeCompetitor(link);
    if (net_id == 827) // Fat fish family
        return generateFatFishCompetitor(link, config.net827_pop_type);
    if (net_id == 1096)
        return generatePoppedCone(link, {config.net1096_removed_face}, pop_description.at(net_id).second);
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

// Tries to find a competitor sufficient to show that all nets in the family
// passed are not minimizing.
// The net is specified via net_id and, if the net is not (almost) rigid, via
// some additional real parameters.
// One can decide to skip running the surface evolver and use the cached
// competitors.
// In any case, this function checks whether the competitors are valid.
// If the function finds a valid competitor, it saves it and returns the (positive)
// difference between the cone and the competitor.
// Otherwise, it returns an error.
ErrorOr<double> excludeFamilyOfNets(int net_id, vector<Interval> net_params,
                                    bool run_surface_evolver,
                                    FindCompetitorConfig config) {
    // Returns a NetEmbeddedWithPins representing the whole family.
    ApproximateGeodesicNetAntipodalPins family = embedNetWithPins(net_id, net_params, false);
    // Chooses a specific net (usually "in the middle of" the family).
    GeodesicNetAntipodalPins link =
            makeNetDeterministic(embedNetWithPins(net_id, net_params, true));
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

// Returns gap and guarantees that the given fat fish net admits a competitor
// (in B_1) which improves the area at least by gap.
// If gap = -1, it means that no gap was found.
double processFatFish(double theta0, double theta1) {
    // These configurations determine the popped cone to use and the
    // surface evolver settings to use.
    // The first parameters (in {-1, 0, 1}) determines the popped cone, while
    // the second determines whether we should execute a "quick and shallow" run
    // of the surface evolver or a "slow and deep" run.
    // We run six different configurations in order and we stop as soon as
    // we find a satisfactory competitor.
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
    cout << " gap = " << area_gap;

    return area_gap;
}

// Consider the family of fat fish corresponding to the interval of parameters
// theta0 and theta1.
//  - theta0 controls the position of the pins (geodesic length of an edge to a pin)
//  - theta1 controls the derived-once net (geodesic length of an edge of the net)
//
// This function returns:
//  - eps > 0 (= 1e-10) if there is no minimizing fat fish in the
//    family (i.e., the parameters are out of bound) but we don't have
//    produced any competitor to prove it.
//  - gap > 0 (a positive number) if one of the nets in the family admits a
//    competitor that improves the area by gap and these competitor is sufficient,
//    employing the stability results, to conclude that no net in the family
//    is minimizing.
//  - -1.0 if the function proves nothing.
double proveInstabilityFatFishRecursive(Interval theta0, Interval theta1) {
    static double area_processed = 0;

    if (theta0.len() + theta1.len() < 1e-4)
        return -1.0;

    double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1.r)));
    // Here we are using the a-priori estimate on the distance of a pin
    // from an edge that is not adjacent to an edge adjacent to the pin itself.
    // This allows us to enforce the -0.01 restriction.
    theta0.intersect(Interval(-1000, upperbound_theta0 - 0.01));
    if (theta0.empty() or theta1.empty())
        return 1e-10;

    cout << theta0 << " " << theta1 << ", " << (100*area_processed / 0.36) << "%, ";

    // If any net in the family admits a competitor with gap > desired_area_gap,
    // then automatically no net in the family can be minimizing.
    // If the net admitting the competitor has parameter corresponding exactly
    // to the midpoints of theta0 and theta1, then desired_area_gap/2 is
    // sufficient.
    double desired_area_gap = 10 * theta0.len() + 4 * theta1.len() + 1e-4;

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

    // We recur on four subrectangles.
    // If the competitor produced while solving any of these subrectangles
    // has gap > desired_area_gap, we can immediately terminate.
    // This observation, as well as many others in processFatFish
    // is necessary to show the non-minimality of all fat fish. Many of these
    // cheap tricks are due to some strange behaviors of the surface evolver
    // in pathological cases of fat fish.
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

void proveInstabilityFatFish() {
    cout << "Processing fat fish." << endl;
    const double sigma = sqrt(2. / 3);
    double upperbound_theta1 = acos(sigma);
    double res = proveInstabilityFatFishRecursive(
            Interval(0, 1000),
            Interval(0, upperbound_theta1-0.01));
    assert(res > 0);
    cout << "OK" << endl << endl;
}


void proveInstabilityRigidNets() {
    cout << "Processing rigid nets." << endl;
    cout << "net_id: area gap" << endl;
    for (int net_id: rigid_nets) {
        cout << net_id << ": \t";
        cout.flush();

        auto maybe_proof = excludeFamilyOfNets(net_id, {}, false, FindCompetitorConfig());
        if (!maybe_proof.isOk())
            cout << maybe_proof.getError() << endl;
        assert(maybe_proof.isOk());
        cout << endl;
    }
    cout << endl;
}

// If this function terminates, it is guaranteed that there is a finite set family
// of subintervals I_1, .., I_k such that:
// - They cover I, that is I \subseteq \bigcup_j I_j
// - It holds f(I_j) = true for all j.
void coverInterval(Interval I, function<bool(Interval)> f) {
    if (f(I))
        return;
    assert(I.len() > 1e-5);
    Interval left = {I.l, I.mid()};
    Interval right = {I.mid(), I.r};
    left.intersect(I);
    right.intersect(I);
    coverInterval(left, f);
    coverInterval(right, f);
};

void proveInstabilityRigidNetsWithFlexiblePins() {
    cout << "Processing rigid nets with flexible pins." << endl;
    for (int net_id: rigid_nets_flexible_pins) {
        cout << net_id << endl;

        ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
        assert(metrized_net_in.is_open());
        ApproximatelyMetrizedNet net = parseApproximatelyMetrizedNet(metrized_net_in);
        // param denotes the length of an edge to a pin.
        Interval param_region = validPinEdgeLength(net);
        if (net_id == 1096) {
            // We use the fact that, under suitable assumptions, the pins
            // cannot be too close to an edge that is not adjacent to an
            // edge adjacent to the pin itself.
            // This condition is non-empty only in the case of net 1096 and
            // it is necessary to succesfully employ the stability results that
            // allow us to discretize the problem.
            param_region.l += 0.005;
            param_region.r -= 0.005;
        }
        clog << "Param region: " << param_region << endl;

        // Returns whether the given interval of parameters contain only
        // non minimal cones.
        auto try_interval = [&](Interval param) {
            if (param.len() > 0.01)
                return false;
            clog << "Processing " << param << endl;
            FindCompetitorConfig config;
            // For net 1096, we need to use two different popped cones depending
            // on param.
            config.net1096_removed_face = (net_id == 1096 and param.mid() > 0.38) ? 8 : 6;
            return excludeFamilyOfNets(net_id, {param}, true, config).isOk();
        };

        coverInterval(param_region, try_interval);

        cout << "OK" << endl << endl;
    }
}

void proveInstabilityGoggles() {
    cout << "Processing goggles." << endl;
    // Here we are using the a-priori lower bound on the length of edges.
    Interval param_region = {0., acos(1.0/3) / 2 - 0.005};
    clog << "Param region: " << param_region << endl;

    // Returns whether the given interval of parameters contain only
        // non minimal cones.
    auto try_interval = [&](Interval param) {
        if (param.len() > 0.01)
            return false;
        clog << "Processing " << param << endl;
        FindCompetitorConfig config;
        config.net1029_cube_family = false;
        // The second parameter Interval(0, 0) is necessary to transmit the
        // information that we are interested in the GOGGLES instead of the
        // cube. Indeed these two families have the same once-derived net (1029).
        return excludeFamilyOfNets(GOGGLES_ID, {param, Interval(0, 0)}, true, config).isOk();
    };

    coverInterval(param_region, try_interval);

    cout << "OK" << endl << endl;
}

// This function can be used to check whether the surface evolver manages
// to produce a better competitor for the cube.
// For 0 <= q < 0.1 it seems to work, for q > 0.1 it does not.
void analyzeCubes(double q) {
    FindCompetitorConfig config;
    config.net1029_cube_family = false;
    vector<Interval> params = {Interval(q), Interval(1, 1)};
    string name = generateEmbeddingFilename(1029, params);
    cout << name << endl;
    cout << excludeFamilyOfNets(GOGGLES_ID, params, true, config) << endl;
}

int main() {
    proveInstabilityRigidNets();
    proveInstabilityRigidNetsWithFlexiblePins();
    proveInstabilityGoggles();
    proveInstabilityFatFish();

    return 0;
}