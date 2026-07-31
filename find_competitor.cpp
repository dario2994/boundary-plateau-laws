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
    vector<Point> embedding = link.deterministicEmbedding();
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
    link.printInfo();
    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedConeCases(net_id, link, config);

    const string& embedding_name = generateEmbeddingFilename(net_id, net_params);
    const string& net_directory = format("data/surface_evolver/{}", net_id);
    const string& filename_prefix = format("{}/{}", net_directory, embedding_name);

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

    auto maybe_validation_error = isCompetitorValid(competitor, link);
    if (!maybe_validation_error.isOk())
        return ErrorMessage(maybe_validation_error.getError());

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

// If this function terminates, it is guaranteed that there is a finite set S of
// points strictly inside (l, r) such that:
//  [l, r] \subseteq \bigcup_{s\in S} interval_covered(s)
// where interval_covered(s) = [s-distance_covered(s), s+distance_covered(s)].
void coverInterval(double l, double r, function<double(double)> distance_covered) {
    // We try to avoid to call distance_covered(s) with values of s very close
    // to l or r (because it can lead to failures in our problem).
    // This is the reason why the code is a bit involved.
    double new_l  = (l + r) / 2;
    while (distance_covered(new_l) < new_l - l)
        new_l = (new_l + l) / 2;

    double new_r  = (l + r) / 2;
    while (distance_covered(new_r) < r - new_r)
        new_r = (new_r + r) / 2;

    clog << "new_l = " << new_l << ", new_r = " << new_r << endl;

    // We have covered [l, new_l] and [new_r, r].
    // Now we focus on the interval [new_l, new_r].
    double param = new_l;
    while (param < new_r) {
        param += distance_covered(param);
    }
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
    // if (ssize(error_messages) == 6) {
    //     for (string msg: error_messages)
    //         cout << msg << endl << endl;
    //     assert(0);
    // }
    cout << " " << area_gap << " (" << ssize(error_messages) << ")";

    return area_gap;
}

double area_processed = 0;
// theta0 controls the position of the pins
// theta1 controls the derived-once net
double proveInstabilityFatFishRecursive(Interval theta0, Interval theta1) {
    if (theta0.len() + theta1.len() < 1e-4)
        return -1.0;

    double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1.r)));
    theta0.intersect(Interval(-1000, upperbound_theta0 - 0.01));
    if (theta0.empty())
        return true;
    cout << theta0 << " " << theta1 << ", " << (100*area_processed / 0.36) << "%, ";

    double desired_area_gap =
        8 * theta0.len() + 2*sqrt(3) * theta1.len() + 1e-3;
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
        cout << net_id << ": \t";
        cout.flush();
        auto maybe_proof = processNet(net_id, {}, true, FindCompetitorConfig());
        if (!maybe_proof.isOk())
            cout << maybe_proof.getError() << endl;
        assert(maybe_proof.isOk());
        double area_gap = maybe_proof.getContent();
        cout << area_gap << endl;
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
        clog << "Param region: " << param_region << endl;

        auto get_area_gap = [&](double param) {
            FindCompetitorConfig config;
            config.net1096_removed_face = (net_id == 1096 and param > 0.38) ? 8 : 6;
            return valueOrAssert(processNet(net_id, {param}, true, config));
        };

        auto distance_covered = [&](double param) {
            double area_gap = get_area_gap(param);
            clog << "area_gap(" << param << ") = " << area_gap << endl;
            return area_gap / 10.;
        };

        coverInterval(param_region.l, param_region.r, distance_covered);
        cout << "OK" << endl << endl;
    }
}

void proveInstabilityGoggles() {
    cout << "Processing goggles." << endl;
    Interval param_region = {0., acos(1.0/3) / 2};
    clog << "Param region: " << param_region << endl;

    auto get_area_gap = [&](double param) {
        FindCompetitorConfig config;
        config.net1029_cube_family = false;
        return valueOrAssert(processNet(GOGGLES_ID, {0., param}, true, config));
    };

    auto distance_covered = [&](double param) {
        double area_gap = get_area_gap(param);
        clog << "area_gap(" << param << ") = " << area_gap << endl;
        return area_gap / 10.;
    };

    coverInterval(param_region.l, param_region.r, distance_covered);

    return;
}

int main() {
    // proveInstabilityRigidNets();
    // proveInstabilityRigidNetsWithFlexiblePins();
    proveInstabilityGoggles();
    // proveInstabilityFatFish();

    return 0;
}

// The kind of estimate we actually need is:

// Sia E una rete geodetica con due pins antipodali.
// Sia \tilde E un'altra rete geodetica con due pins antipodali.

// Creiamo un competitor valido per \tilde E con un'area gap di g.

// Se E è una rete valida, allora g < dist(E, \tilde E)

//
// Rigid case
//
// Goal: Any net compatible with the approximate embedding Q is not stable.
// Method:
// Take the center of each point in the approximate embedding.
// In this way we obtain a geodesic net N.
// We find a competitor for N with area gap g.
// Now, assume that there is a stable net in Q, let it be \tilde N.
// If dist(N, \tilde N) is small compared to g, then this is impossible.


// Lemma: Let $E$ be a cone containing the line l whose link is a geodesic net with two pins.
//        Assume that E admits a competitor with Area(F(E)) <= Area(E) - \eps.
//        Let $\tilde E$ be another cone whose link is a geodesic net with two pins.
//        There is a competitor for $\tilde E$ such that
//        Area(\tilde F(\tilde E)) <= Area(\tilde E) - \eps + F(dist(E, \tilde E).

//
// Lemma: Let $E$ be a cone containing the line l whose link is a geodesic net with two pins.
//        Assume that E admits a competitor with Area(competitor) <= Area(E) - \eps.
//        Let $\tilde E$ be another cone whose link is a geodesic net with two pins.
//
//
//
//
//
//
//
