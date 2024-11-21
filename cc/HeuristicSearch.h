#pragma once
#include <NTL/LLL.h>
#include <fplll/fplll.h>

namespace ntl {

double lll(NTL::Vec<NTL::ZZ> &, int, int);
double hkz(NTL::Vec<NTL::ZZ> &, int, int);
} // namespace ntl

struct longvec {
    NTL::Vec<NTL::ZZ> v;
    double l;
};

class HeuristicSearch {
    using svp = double (*)(NTL::Vec<NTL::ZZ> &, int, int);
    int N;
    svp F;

  public:
    HeuristicSearch(int N, std::string s);
    HeuristicSearch(int N, svp F);

    longvec cube(int p, double r);
    longvec simplex(int p, double r);
    longvec diagonal(int p, double r1, double r2);
};

namespace fplll {
typedef fplll::ZZ_mat<mpz_t> mat_zz;
typedef fplll::Z_NR<mpz_t> zz;

double lll(NTL::Vec<NTL::ZZ> &, int, int);
double hkz(NTL::Vec<NTL::ZZ> &, int, int);

}; // namespace fplll
