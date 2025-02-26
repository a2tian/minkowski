#pragma once
#include <NTL/LLL.h>

#include <vector>

namespace ntl {
NTL::RR lll(NTL::Vec<NTL::ZZ> &, NTL::ZZ p, long n);
NTL::RR hkz(NTL::Vec<NTL::ZZ> &, NTL::ZZ p, long n);
}  // namespace ntl

struct longvec {
    NTL::Vec<NTL::ZZ> v;
    NTL::RR l;
};

class HeuristicSearch {
    using svp = NTL::RR (*)(NTL::Vec<NTL::ZZ> &, NTL::ZZ p, long n);
    svp F;

   public:
    HeuristicSearch(svp F);
    HeuristicSearch(std::string s = "hkz");

    longvec one(long n, NTL::ZZ p, NTL::Vec<NTL::ZZ> a);
    longvec cube(long n, NTL::ZZ p, double r = 0.25);
    longvec simplex(long n, NTL::ZZ p, double r = 0.25);
    longvec diagonal(long n, NTL::ZZ p, double r1 = 0.25, double r2 = 0.25);

    // TODO: Implement in NTL
    // longvec center(long n, NTL::ZZ p, NTL::ZZ w, std::vector<NTL::ZZ> c, double r1 = 0.25, double r2 = 0.25);
};
