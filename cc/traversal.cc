#include "traversal.h"

#include <gmp.h>

using namespace NTL;
using namespace std;

/**
 * @brief Compute q given p
 */
double q(int p, int N) {
    return pow(p, (N - 1.) / N);
}

/**
 * @brief Compute the Minkowski lower bound for a shortest lattice vector given N
 */
double minkowski(int N) {
    return std::pow(2 * std::tgamma(N / 2. + 1) / std::pow(M_PI, N / 2.), 1. / N);
}

namespace traversal {
/**
 * @brief Increment a vector in U(p)
 */
void cube(Vec<ZZ> &a, int b, int n) {
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
void simplex(Vec<ZZ> &a, int, int n) {
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
void diagonal(Vec<ZZ> &a, double r, int n) {
    for (int i = n - 1; i > 1; i--) {
        if (a[i] < a[i - 1]) {
            a[i]++;
            return;
        }
        a[i] = (int)((1 - r) * conv<int>(a[i - 1]));
    }
    a[1]++;
}

void cdiagonal(Vec<ZZ> &a, double r, Vec<int> b1, Vec<int> b2, int n) {
    for (int i = n - 1; i > 1; i--) {
        if (a[i] < min(conv<int>(a[i - 1]), b2[i - 1])) {
            a[i]++;
            return;
        }
        a[i] = max((int)((1 - r) * conv<int>(a[i - 1])), b1[i - 1]);
    }
    a[1]++;
}
} // namespace traversal