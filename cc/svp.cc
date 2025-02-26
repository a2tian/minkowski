#include "HeuristicSearch.h"
#include "vector.h"

using namespace NTL;
using namespace std;

namespace ntl {
    
using svp = RR (*)(Mat<ZZ> &, long);

RR EXACT_LLL(Mat<ZZ> &B, long) {
    ZZ det2;
    LLL(det2, B, 99, 100);
    return norm(B[0]);
}

RR HKZ(Mat<ZZ> &B, long n) {
    BKZ_FP(B, 0.99, n);
    return norm(B[0]);
}

/**
 * @brief Compute the length of the shortest vector
 * in the R-reduced basis of B(p, a)
 *
 * @param p a non-negative integer
 * @param a a vector of type NTL::Vec<NTL::ZZ>
 * @param R a pointer to a function that performs a lattice basis reduction
 * @param N the dimension of the lattice
 */
RR shortest_vector(Vec<ZZ> &a, svp f, ZZ p, long n) {
    Mat<ZZ> B = ident_mat_ZZ(n) * p;
    B[0] = a;
    return f(B, n);
}

/**
 * @brief Compute the length of the shortest vector in the
 * LLL-reduced basis of B(p, a)
 *
 * @param a a vector of type NTL::Vec<NTL::ZZ>
 * @param p a non-negative integer
 * @param N the dimension of the lattice
 */
RR lll(Vec<ZZ> &a, ZZ p, long n) {
    return shortest_vector(a, EXACT_LLL, p, n);
}

/**
 * @brief Compute the length of the shortest vector in the
 * KZ-reduced basis of B(p, a)
 *
 * @param a a vector of type NTL::Vec<NTL::ZZ>
 * @param p a non-negative integer
 * @param N the dimension of the lattice
 */
RR hkz(Vec<ZZ> &a, ZZ p, long n) {
    return shortest_vector(a, HKZ, p, n);
}
} // namespace ntl
