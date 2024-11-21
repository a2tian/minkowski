#include "vector.h"

using namespace NTL;
using namespace std;

namespace ntl {
/**
 * @brief Compute the norm of a vector
 */
double norm(Vec<ZZ> &v) {
    double res = 0;
    for (int i = 0; i < v.length(); i++) {
        res += conv<int>(v[i] * v[i]);
    }
    return sqrt(res);
}
} // namespace ntl

namespace fplll {
double mat_norm(mat_zz &M, int n) {
    int res = 0;

    for (int i = 0; i < n; i++) {
        int el = M[0][i].get_si();
        res += el * el;
    }
    return std::sqrt(res);
}

/**
 * @brief Generates pI, where I is the nxn identity matrix.
 *
 * @param M the matrix to be modified
 * @param n dimension
 * @param p (default 1) scalar
 */
void gen_identity(mat_zz &M, int n, int p = 1) {
    M.resize(n, n);
    for (int i = 0; i < n; i++) {
        M(i, i) = p;
    }
}

/**
 * @brief Initializes a matrix for enumeration of U(p)
 *
 * @param M the matrix to be modified
 * @param p scalar
 * @param n dimension
 * @param i (default 0) the value to be assigned to the second element of the first row
 */
void init_matrix(mat_zz &M, Vec<ZZ> &a, int p, int n) {
    M.resize(n, n);
    for (int i = 0; i < n; i++) {
        M(i, i) = p;
        M(0, i) = conv<int>(a[i]);
    }
}

/**
 * @brief Initializes the lattice B(p, a) and increments a
 *
 * @param M the matrix to be modified
 * @param p scalar
 * @param n dimension
 * @param i (default 0) the value to be assigned to the second element of the first row
 */
void init_and_increment(mat_zz &M, std::vector<int> &a, int p, int b, int n) {
    M.resize(n, n);
    int flag = 0;
    for (int i = n - 1; i > 0; i--) {
        M(i, i) = p;
        M(0, i) = a[i];
        a[i] = (a[i] < b - 1) ? a[i] + (++flag) : flag * a[i];
    }
    a[1] += (flag == 0) ? 1 : 0;
}
} // namespace fplll