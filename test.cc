#include <algorithm>
#include <chrono>
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <random>

#include "HeuristicSearch.h"
#include "traversal.h"

using namespace std;
using namespace NTL;
namespace fs = std::experimental::filesystem;

void gen_random_vector(Vec<ZZ> &a, ZZ p, long n) {
    a.SetLength(n - 1);
    for (int i = 0; i < n - 1; i++) {
        a[i] = RandomBnd(p);
    }
}

void gen_random_vector_center(Vec<ZZ> &a, ZZ p, long n, RR frac) {
    a.SetLength(n - 1);
    ZZ r = conv<ZZ>(conv<RR>(p) * frac);

    for (int i = 0; i < n - 1; i++) {
        a[i] = RandomBnd(p) - r;
    }
}

void make_center(Vec<ZZ> &res, Vec<ZZ> &a, long n, ZZ p, RR frac) {
    res.SetLength(n - 1);
    ZZ r = conv<ZZ>(conv<RR>(p) * frac);

    for (int i = 0; i < n - 1; i++) {
        res[i] = r + a[i];
    }
}

void random_m_decreasing(long n, int m, int rand_offset, int init_offset, ofstream &out) {
    Vec<ZZ> a;
    ZZ r = power2_ZZ(n + rand_offset), p = NextPrime(power2_ZZ(n + init_offset));
    gen_random_vector(a, r, n);

    HeuristicSearch hs;
    for (int i = 0; i < m; i++) {
        out << p << " " << hs.one(n, p, a).l / q(p, n) << endl;
        p = NextPrime(p * 2);
    }
}

void random_m_decreasing(long n, int m, int rand_offset, int init_offset, RR frac, ofstream &out) {
    Vec<ZZ> a;
    ZZ r = power2_ZZ(n + rand_offset), p = NextPrime(power2_ZZ(n + init_offset));
    gen_random_vector_center(a, p, n, frac);

    HeuristicSearch hs;
    for (int i = 0; i < m; i++) {
        Vec<ZZ> res;
        make_center(res, a, n, p, frac);
        out << p << " " << hs.one(n, p, res).l / q(p, n) << endl;
        p = NextPrime(p * 2);
    }
}

// void offset(Vec<ZZ> &res, Vec<ZZ> &a, ZZ p, double r, long n) {
//     res.SetLength(n - 1);
//     for (int i = 0; i < n - 1; i++) {
//         res[i] = a[i] + p * r;
//     }
// }

void clamp0p(ZZ &val, ZZ &p) {
    if (val < 0) val = 0;
    if (val >= p) val = p - 1;
}

// std::vector<double> spsaOptimize(
//     long n,
//     int maxIters,
//     double a,      // step-size coefficient
//     double c,      // perturbation-size coefficient
//     double alpha,  // exponent for step-size decay
//     double gamma,  // exponent for perturbation decay
//     double A       // stability constant in denominator of step-size
// ) {
//     // Initialize x.
//     // We'll set x[0] = 1.0 explicitly and randomly init the rest in [0,1].

//     ZZ p = NextPrime(power2_ZZ(n));
//     Vec<ZZ> x;
//     gen_random_vector(x, p, n);

//     HeuristicSearch hs;

//     // Evaluate initial value
//     RR bestVal = hs.one(n, p, x).l;
//     Vec<ZZ> bestX = x;

//     // Set up random sign generator for SPSA:
//     static std::mt19937 rng(54321);
//     std::uniform_int_distribution<int> signDist(0, 1);
//     // signDist(rng) gives 0 or 1; we'll map that to +1 or -1.

//     // Main SPSA loop
//     for (int k = 0; k < maxIters; ++k) {
//         // Compute step sizes for this iteration
//         double a_k = a / std::pow(k + 1.0 + A, alpha);
//         double c_k = c / std::pow(k + 1.0, gamma);

//         // Construct random sign vector Delta in dimension n-1 (x[0] is fixed)
//         std::vector<double> Delta(n, 0.0);
//         for (int i = 1; i < n; ++i) {
//             int s = signDist(rng);
//             Delta[i] = (s == 0) ? +1.0 : -1.0;
//         }

//         // Form x_plus and x_minus
//         Vec<ZZ> x_plus = x;
//         Vec<ZZ> x_minus = x;
//         for (int i = 1; i < n; ++i) {
//             x_plus[i] = clamp01(x[i] + c_k * Delta[i]);
//             x_minus[i] = clamp01(x[i] - c_k * Delta[i]);
//         }

//         // Evaluate sphere packing "density" at x_plus and x_minus
//         RR f_plus = hs.one(n, p, x_plus).l;
//         RR f_minus = hs.one(n, p, x_minus).l;

//         // Compute approximate gradient:
//         //   (f_plus - f_minus)/(2*c_k) * Delta
//         std::vector<double> grad(n, 0.0);
//         double diff = (f_plus - f_minus) / (2.0 * c_k);
//         for (int i = 1; i < n; ++i) {
//             grad[i] = diff * Delta[i];
//         }

//         // Update x with step size a_k
//         std::vector<double> x_new = x;
//         for (int i = 1; i < n; ++i) {
//             x_new[i] = clamp01(x[i] + a_k * grad[i]);
//         }

//         // Evaluate new point
//         double f_new = spherePackingDensity(n, x_new);

//         // Keep if better
//         if (f_new > bestVal) {
//             bestVal = f_new;
//             bestX = x_new;
//         }

//         // Update current x for next iteration
//         x = x_new;
//     }

//     // Return the best found solution
//     return bestX;
// }

int main() {
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    int m = 3;

#pragma omp parallel for
    for (int n = 53; n <= 64; n++) {
        for (int rand_offset = 0; rand_offset < 1; rand_offset++) {
            string dir = "out/corner/" + to_string(n);
            if (!fs::exists(dir)) {
                fs::create_directories(dir);
            }
            string file = dir + "/" + to_string(rand_offset) + "_fast.txt";
            ofstream out(file);
            auto t1 = high_resolution_clock::now();
            random_m_decreasing(n, m, rand_offset, rand_offset, out);
            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<milliseconds>(t2 - t1);

            out << ms_int.count() << "ms" << endl;
            out.close();
        }

        // for (int rand_offset = -n / 2; rand_offset <= 2 * n; rand_offset += 2) {
        //     string dir = "out/1-2/" + to_string(n);
        //     if (!fs::exists(dir)) {
        //         fs::create_directories(dir);
        //     }
        //     string file = dir + "/" + to_string(rand_offset) + ".txt";
        //     ofstream out(file);
        //     auto t1 = high_resolution_clock::now();
        //     random_m_decreasing(n, m, rand_offset, rand_offset, RR(0.5), out);
        //     auto t2 = high_resolution_clock::now();
        //     auto ms_int = duration_cast<milliseconds>(t2 - t1);

        //     out << ms_int.count() << "ms" << endl;
        //     out.close();
        // }
    }
}