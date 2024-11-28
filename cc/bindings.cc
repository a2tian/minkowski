#include "HeuristicSearch.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_minkowski, m) {
    m.doc() = "A module for finding dense lattice packings of hyperspheres";

    nb::bind_vector<std::vector<int>>(m, "intvec");

    nb::class_<NTL::Vec<NTL::ZZ>>(m, "Vec_ZZ")
        .def(nb::init<>())
        .def("__init__", [](std::vector<int> v) {
            NTL::Vec<NTL::ZZ> vec;
            for (auto i : v)
                vec.append(NTL::ZZ(i));
            return vec;
        })
        .def("__getitem__", [](const NTL::Vec<NTL::ZZ> &v, size_t i) {
            if (i >= v.length()) throw std::out_of_range("Index out of range");
            return NTL::conv<int>(v[i]); // Convert ZZ to Python int
        })
        .def("__repr__", [](const NTL::Vec<NTL::ZZ> &v) {
            std::ostringstream oss;
            oss << "Vec(" << v << ")";
            return oss.str();
        })
        .def("__len__", [](const NTL::Vec<NTL::ZZ> &v) { return v.length(); });

    nb::class_<longvec>(m, "longvec")
        .def(nb::init<>())
        .def_ro("l", &longvec::l)
        .def_ro("v", &longvec::v)
        .def("__repr__", [](const longvec &lv) {
            std::ostringstream oss;
            oss << "shortest vector: " << lv.v << "; length of shortest vector: " << lv.l;
            return oss.str();
        });

    nb::class_<HeuristicSearch>(m, "HeuristicSearch")
        .def(nb::init<std::string>(), "svp"_a = "hkz")
        .def("one", &HeuristicSearch::one, "n"_a, "p"_a, "a"_a, "Find the length of the shortest vector in the lattice B(p, a)")
        .def("cube", &HeuristicSearch::cube, "n"_a, "p"_a, "r"_a = 0.25, "Search U(p) in a cube of width p * r")
        .def("simplex", &HeuristicSearch::simplex, "n"_a, "p"_a, "r"_a = 0.25, "Search U(p) in the simplex of width p*r")
        .def("diagonal", &HeuristicSearch::diagonal, "n"_a, "p"_a, "r1"_a = 0.25, "r2"_a = 0.25, "Search U(p) in the diagonal of width p*r1 and radius r2")
        .def("center", &HeuristicSearch::center, "n"_a, "p"_a, "w"_a, "c"_a, "r1"_a = 0.25, "r2"_a = 0.25);
}