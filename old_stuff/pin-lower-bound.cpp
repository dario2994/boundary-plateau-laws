#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>

// Basic 3D Vector structure
struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

double norm(Vec3 v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double area(Vec3 p1, Vec3 p2, Vec3 p3) {
    return 0.5 * norm(cross(p2 - p1, p3 - p1));
}

int main() {
    // Fixed Parameters
    const double L = 0.2;
    const double r = 0.8;
    const Vec3 O = {0.0, 0.0, 0.0};
    const Vec3 B = {L, 0.0, 1.0};
    const Vec3 C = {-L, 0.0, 1.0};
    
    // Grid Parameters
    const double delta_min = 0.0;
    const double delta_max = 0.011;
    const double delta_step = 0.0001;
    
    const double alpha_min = M_PI / 4.0; 
    const double alpha_max = M_PI / 2.0;
    const double alpha_step = 0.001;

    // Lipschitz Constants and Safety Margin
    const double K_delta = 4.0;
    const double K_alpha = 1.3;
    const double max_variation = (K_delta * (delta_step / 2.0)) + (K_alpha * (alpha_step / 2.0));

    double global_min_F = 1e9;
    double min_delta = 0.0;
    double min_alpha = 0.0;

    std::cout << "Starting grid evaluation..." << std::endl;
    std::cout << "Alpha range: [PI/4, PI/2], Delta range: [0, 0.011]" << std::endl;

    for (double alpha = alpha_min; alpha <= alpha_max + 1e-9; alpha += alpha_step) {
        // Calculate stationary points for this alpha
        Vec3 A = {L * std::cos(alpha), L * std::sin(alpha), 1.0};
        
        // Explicit O(1) evaluation of the Fermat-Torricelli point D(alpha)
        double xD = L * (0.5 + (1.0 / std::sqrt(3.0)) * std::sin(M_PI / 3.0 - alpha));
        double yD = L * (-1.0 / (2.0 * std::sqrt(3.0)) + (1.0 / std::sqrt(3.0)) * std::cos(M_PI / 3.0 - alpha));
        Vec3 D = {xD, yD, 1.0};
        
        Vec3 Ap = A * r;
        Vec3 Bp = B * r;
        Vec3 Dp = D * r;

        for (double delta = delta_min; delta <= delta_max + 1e-9; delta += delta_step) {
            Vec3 M = {delta * std::cos(alpha), delta * std::sin(alpha), 1.0};
            Vec3 Mp = M * r;

            // Base Area
            double G = area(B, C, O) + area(A, M, O);

            // Competitor Area (using the explicitly defined 10-triangle complex)
            double H = area(B, C, M) + area(M, C, O) +
                       area(Ap, Dp, O) + area(Bp, Dp, O) + area(Mp, Dp, O) +
                       area(A, Ap, Dp) + area(A, M, Dp) + area(M, Mp, Dp) +
                       area(B, M, Dp) + area(B, Bp, Dp);

            double F = G - H;

            if (F < global_min_F) {
                global_min_F = F;
                min_delta = delta;
                min_alpha = alpha;
            }
        }
    }

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n--- Grid Evaluation Complete ---" << std::endl;
    std::cout << "Minimum F found on grid: " << global_min_F << std::endl;
    std::cout << "Occurs at delta = " << min_delta << ", alpha = " << min_alpha << " radians" << std::endl;
    
    std::cout << "\n--- Rigorous Gap Analysis ---" << std::endl;
    std::cout << "Maximum Lipschitz variation within cell: " << max_variation << std::endl;
    std::cout << "Guaranteed minimum bound (Grid Min - Variation): " << global_min_F - max_variation << std::endl;

    if (global_min_F - max_variation > 0) {
        std::cout << "\nSUCCESS: The continuous function F is rigorously proven to be strictly positive." << std::endl;
    } else {
        std::cout << "\nFAILURE: The grid resolution is too coarse, or F crosses zero." << std::endl;
    }

    return 0;
}