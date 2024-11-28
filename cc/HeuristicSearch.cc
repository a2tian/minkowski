#include "HeuristicSearch.h"
#include "traversal.h"
#include "vector.h"

#include <omp.h>

using namespace NTL;
using namespace std;

/**
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * HeuristicSearch implementation based on the NTL library.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

struct longvec &cmpvec(longvec &a, longvec &b) {
    return a.l > b.l ? a : b;
}

#pragma omp declare reduction(maxvec:longvec : omp_out = cmpvec(omp_out, omp_in)) initializer(omp_priv = {Vec<ZZ>(), -1})

void init_vector(Vec<ZZ> &a, int i, int n) {
    a.SetLength(n);
    a[0] = 1, a[1] = i;
}

HeuristicSearch::HeuristicSearch(svp F) : F(F) {}

HeuristicSearch::HeuristicSearch(string s) {
    if (s == "lll")
        F = ntl::lll;
    else if (s == "hkz")
        F = ntl::hkz;
    else
        throw invalid_argument("Invalid argument: " + s + ". Must be 'lll' or 'hkz'.");
}

/**
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Note: the following functions are non-parallelized versions.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Vec<ZZ> HeuristicSearch::cube(int p, double r) {
//     int n = N, b = p * r;
//     svp f = F;
//     struct longvec lv = {Vec<ZZ>(), -1};

//     for (int i = 0; i < b; i++) { // iterate over first dimension
//         Vec<ZZ> a;
//         init_vector(a, i, n);
//         while (a[1] == i) {
//             double l = f(a, p, n);
//             if (l > lv.l) lv = {a, l};
//             traversal::cube(a, b, n);
//         }
//     }
//     return lv.v;
// }

// Vec<ZZ> HeuristicSearch::simplex(int p, double r) {
//     int n = N, b = p * r;
//     svp f = F;

//     struct longvec lv = {Vec<ZZ>(), -1};

//     Vec<ZZ> a;
//     init_vector(a, 0, n);
//     while (a[1] < b) {
//         traversal::simplex(a, b, n);
//         double l = f(a, p, n);
//         if (l > lv.l) lv = {a, l};
//     }
//     return lv.v;
// }

// Vec<ZZ> HeuristicSearch::diagonal(int p, double r1, double r2) {
//     int n = N, b = p * r1;
//     svp f = F;

//     struct longvec lv = {Vec<ZZ>(), -1};

//     Vec<ZZ> a;
//     init_vector(a, 0, n);

//     while (a[1] < b) {
//         traversal::diagonal(a, r2, n);
//         double l = f(a, p, n);
//         if (l > lv.l) lv = {a, l};
//     }
//     return lv.v;
// }

/**
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * Note: parallelized versions of the following functions via
 *       omp reduction are commented out.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

longvec HeuristicSearch::one(int n, int p, vector<int> a) {
    Vec<ZZ> v;
    init_vector(v, 0, n);
    for (int i = 0; i < n - 1; i++)
        v[i + 1] = a[i];
    return {v, F(v, p, n)};
}

longvec HeuristicSearch::cube(int n, int p, double r) {
    svp f = F;
    int b = p * r;
    struct longvec lv = {Vec<ZZ>(), -1};
#pragma omp parallel for reduction(maxvec : lv)
    for (int i = 0; i < b; i++) { // iterate over first dimension
        Vec<ZZ> a;
        init_vector(a, i, n);
        while (a[1] == i) {
            double l = f(a, p, n);
            if (l > lv.l) lv = {a, l};
            traversal::cube(a, b, n);
        }
    }
    return lv;
}

longvec HeuristicSearch::simplex(int n, int p, double r) {
    svp f = F;
    Vec<ZZ> a;
    init_vector(a, 0, n);
    int b = p * r;
    char flag = 0;
    struct longvec lv = {Vec<ZZ>(), -1};
#pragma omp parallel shared(flag) reduction(maxvec : lv)
    {
        while (!flag) {
            Vec<ZZ> ap = a;
#pragma omp critical
            {
                traversal::simplex(a, b, n);
                if (a[1] >= b) flag++;
            }
            double l = f(ap, p, n);
            if (l > lv.l) lv = {ap, l};
        }
    }
    return lv;
}

longvec HeuristicSearch::diagonal(int n, int p, double r1, double r2) {
    int b = p * r1;
    svp f = F;
    Vec<ZZ> a;
    init_vector(a, 0, n);
    char flag = 0;
    struct longvec lv = {Vec<ZZ>(), -1};
#pragma omp parallel shared(flag) reduction(maxvec : lv)
    {
        while (!flag) {
            Vec<ZZ> ap = a;
#pragma omp critical
            {
                traversal::diagonal(a, r2, n);
                if (a[1] >= b) flag++;
            }
            double l = f(ap, p, n);
            if (l > lv.l) lv = {ap, l};
        }
    }
    return lv;
}

longvec HeuristicSearch::center(int n, int p, int w, vector<int> c, double r1, double r2) {
    int b = min((int)(p * r1), c[0] + w);
    svp f = F;
    Vec<ZZ> a;
    Vec<int> b1, b2;
    init_vector(a, 0, n);
    b1.SetLength(n - 1), b2.SetLength(n - 1);
    for (int i = 0; i < n - 1; i++) {
        b1[i] = c[i] - w;
        b2[i] = c[i] + w;
        a[i + 1] = max(0, b1[i]);
    }
    char flag = 0;
    struct longvec lv = {Vec<ZZ>(), -1};
#pragma omp parallel shared(flag) reduction(maxvec : lv)
    {
        while (!flag) {
            Vec<ZZ> ap = a;
#pragma omp critical
            {
                traversal::cdiagonal(a, r2, b1, b2, n);
                if (a[1] >= b) flag++;
            }
            double l = f(ap, p, n);
            if (l > lv.l) lv = {ap, l};
        }
    }
    return lv;
}
