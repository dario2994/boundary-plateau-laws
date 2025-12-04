#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"
#include "simplicial_complexes.cpp"
#include "generate_surface_evolver_input.cpp"

// TODO: This should take a competitor instead!
bool isCompetitorValid(ImmersedSimplicialComplex X, NetEmbeddedWithPins link) {
    auto pin1 = link.embedding[link.pins.first];
    auto pin2 = link.embedding[link.pins.second];
    auto maybe_line = X.containsSegment(pin1, pin2);
    if (!maybe_line)
        return false;
    auto line = *maybe_line;

    if (!X.isCollapsible(line))
        return false;

    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;

    // }
    // for (auto [u, v]: link.boundary_edges) {
        auto U = link.embedding[d.from];
        auto V = link.embedding[d.to];
        if (!X.containsSegment(U, V))
            return false;
    }
    return true;
}

void processNet(int net_id, bool run_surface_evolver, vector<double> params={}) {
    NetEmbeddedWithPins link = embedNetWithPins(net_id, params);

    if (!pop_description.contains(net_id)) {
        printInfoForHumanPopping(link);
        return;
    }

    SurfaceEvolverData cone = generateCone(link);
    SurfaceEvolverData popped_cone = generatePoppedCone(link,
        pop_description.at(net_id).first,
        pop_description.at(net_id).second);

    string embedding_name = generateEmbeddingFilename(net_id, params);

    if (run_surface_evolver) {
        assert(0); // TODO
        string filename_prefix = "data/surface_evolver/" + embedding_name;

        saveSurfaceEvolverFormat(
            embedding_name, cone, filename_prefix + "_cone.fe", false);
        saveSurfaceEvolverFormat(
            embedding_name, popped_cone, filename_prefix + "_popped.fe", true);
        // Run the surface evolver and the python script to convert to a reasonable format.
    }

    ifstream competitor_in("data/surface_evolver/" + embedding_name + "_competitor.txt");
    ImmersedSimplicialComplex competitor = parseImmersedSimplicialComplex(competitor_in);

    cout << "#vertices" << endl;
    cout << "Cone: " << cone.X.V() << ", "
         << "Popped cone: " << popped_cone.X.V() << ", "
         << "Competitor: " << competitor.V() << ", ";
    cout << endl;

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


    // pair<int,int> pins = {22, 23};
    // vector<pair<int,int>> boundary_edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 4}, {1, 5}, {5, 6}, {7, 19}, {7, 8}, {8, 9}, {6, 9}, {7, 10}, {10, 11}, {11, 12}, {8, 12}, {13, 20}, {9, 13}, {5, 14}, {0, 14}, {2, 18}, {10, 15}, {13, 14}, {16, 21}, {11, 16}, {3, 17}, {16, 17}, {15, 17}, {15, 18}, {6, 19}, {12, 20}, {4, 21}, {21, 22}, {20, 22}, {19, 23}, {18, 23}, };
    // vector<array<double, 3>> points = {{1, 0, 0}, {0.160836, 0.986981, -5.55112e-17}, {0.0248617, 0.973201, -0.22861}, {0.416452, 0.354083, -0.837373}, {0.96453, -0.131987, -0.228608}, {-0.353093, 0.550021, 0.756837}, {-0.586183, 0.453669, 0.671248}, {-0.955015, 0.140076, -0.26139}, {-0.735364, -0.625231, 0.261392}, {-0.586048, -0.498279, 0.638957}, {-0.85287, 0.125095, -0.506916}, {-0.013706, -0.861885, -0.506918}, {-0.0153461, -0.96511, -0.261393}, {0.353482, -0.65152, 0.671246}, {0.48607, -0.43696, 0.756837}, {-0.723245, 0.337136, -0.602707}, {0.216425, -0.768049, -0.602711}, {0.113592, 0.096581, -0.988822}, {-0.586592, 0.791318, -0.172414}, {-0.761386, 0.642701, 0.0850065}, {0.511874, -0.854845, 0.0849996}, {0.686676, -0.706222, -0.172413}, {0.647757, -0.761847, -7.34244e-06}, {-0.647757, 0.761847, 7.34244e-06}, };

    // ifstream foo_in("data/surface_evolver/182_competitor.txt");
    // ImmersedSimplicialComplex competitor = parseImmersedSimplicialComplex(foo_in);

    // cout << "Area before: " << X.area() << endl;
    // cout << "#vertices = " << competitor.V() << endl;

    // NetEmbeddedWithPins boundary{
    //     .A = constructOneLoop(),
    //     .embedding = points,
    //     .pins = pins,
        // .boundary_edges = boundary_edges
    // };

    // SurfaceEvolverData competitor {
    //     .X = X,
    //     .pins = pins,
    //     .boundary_edges = boundary_edges
    // };
    // NetEmbeddedWithPins link = embedNetWithPins(182, {});
    // cout << isCompetitorValid(competitor, link) << endl;

    // cout << "Area after: " << competitor.area() << endl;
    processNet(182, false, {});
}
