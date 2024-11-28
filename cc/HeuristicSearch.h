#pragma once
#include <NTL/LLL.h>
#include <nanobind/ndarray.h>
#include <vector>

namespace nb = nanobind;

namespace ntl {
double lll(NTL::Vec<NTL::ZZ> &, int p, int n);
double hkz(NTL::Vec<NTL::ZZ> &, int p, int n);
} // namespace ntl

struct longvec {
    NTL::Vec<NTL::ZZ> v;
    double l;
};

class HeuristicSearch {
    using svp = double (*)(NTL::Vec<NTL::ZZ> &, int p, int n);
    using ndarray = nb::ndarray<int, nb::shape<-1>>;
    svp F;

  public:
    HeuristicSearch(svp F);
    HeuristicSearch(std::string s = "hkz");

    longvec one(int n, int p, std::vector<int> a);
    longvec cube(int n, int p, double r = 0.25);
    longvec simplex(int n, int p, double r = 0.25);
    longvec diagonal(int n, int p, double r1 = 0.25, double r2 = 0.25);
    longvec center(int n, int p, int w, std::vector<int> c, double r1 = 0.25, double r2 = 0.25);
};
