#include "traversal.h"

#include <gmp.h>

using namespace NTL;
using namespace std;

/**
 * @brief Compute q given p
 */
RR q(ZZ p, long n) {
    return exp(conv<RR>((n - 1.) / n) * log(conv<RR>(p)));
}

/**
 * @brief Compute the Minkowski lower bound for a shortest lattice vector given N
 */
RR minkowski(long n) {
    return exp(conv<RR>(1. / n) * (log(2) + std::lgamma(n / 2. + 1) - (n / 2.) * log(M_PI)));
}

namespace traversal {
/**
 * @brief Increment a vector in U(p)
 */
void cube(Vec<ZZ> &a, ZZ b, long n) {
    for (int i = n - 1; i > 1; i--) {
        if (a[i] < b - 1) {
            a[i]++;
            return;
        }
        a[i] = 0;
    }
    a[1]++;
}

/**
 * @brief Increment a vector in the symmetric subdivision of U(p)
 */
void simplex(Vec<ZZ> &a, ZZ, long n) {
    for (int i = n - 1; i > 1; i--) {
        if (a[i] < a[i - 1]) {
            a[i]++;
            return;
        }
        a[i] = 0;
    }
    a[1]++;
}

/**
 * @brief Increment a vector in the symmetric subdivision of U(p) within a radius along the diagonal
 */
void diagonal(Vec<ZZ> &a, double r, long n) {
    for (int i = n - 1; i > 1; i--) {
        if (a[i] < a[i - 1]) {
            a[i]++;
            return;
        }
        a[i] = (int)((1 - r) * conv<int>(a[i - 1]));
    }
    a[1]++;
}
}  // namespace traversal

// TODO: Implement in NTL
// void cdiagonal(Vec<ZZ> &a, double r, Vec<ZZ> b1, Vec<ZZ> b2, long n) {
//     for (int i = n - 1; i > 1; i--) {
//         if (a[i] < min(conv<int>(a[i - 1]), b2[i - 1])) {
//             a[i]++;
//             return;
//         }
//         a[i] = max((int)((1 - r) * conv<int>(a[i - 1])), b1[i - 1]);
//     }
//     a[1]++;
// }
// } // namespace traversal