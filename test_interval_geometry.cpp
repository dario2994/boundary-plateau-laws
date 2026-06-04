#include "utils.cpp"
#include "interval_geometry.cpp"


double randomValue(double a, double b) {
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> distrib(a, b);

    return distrib(gen);
}

Interval randomInterval(double a, double b) {
    double l = randomValue(a, b);
    double r = randomValue(a, b);
    if (l > r)
        swap(l, r);
    return Interval{l, r};
}

int main () {
    for (int _ = 0; _ < 1'000'000; _++) {
        Interval X = randomInterval(-2, 2);
        Interval Y = randomInterval(-2, 2);
        Interval domain = randomInterval(0, 2*M_PI-1e-6);
        double theta = randomValue(domain.l, domain.r);
        if (X.l <= cos(theta) and cos(theta) <= X.r and Y.l <= sin(theta) and sin(theta) <= Y.r) {
            Interval I = intersectWithCircle(X, Y, domain);
            assert(I.l <= theta and theta <= I.r);
            assert(isSubset(I, domain));
        }
    }
}
