#include "HeuristicSearch.h"
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <sstream>

namespace nb = nanobind;

NB_MODULE(_minkowski, m) {
    m.doc() = "A module for finding dense lattice packings of hyperspheres";
    nb::class_<NTL::Vec<NTL::ZZ>>(m, "VecZZ")
        .def(nb::init<>())
        .def("__getitem__", [](const NTL::Vec<NTL::ZZ> &v, size_t i) {
            if (i >= v.length()) throw std::out_of_range("Index out of range");
            return NTL::conv<int>(v[i]); // Convert ZZ to Python int
        })
        .def("__repr__", [](const NTL::Vec<NTL::ZZ> &v) {
            std::ostringstream oss;
            oss << "Vec_ZZ(" << v << ")";
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
        .def(nb::init<int, std::string>())
        .def("cube", &HeuristicSearch::cube, nb::call_guard<nb::gil_scoped_acquire>())
        .def("simplex", &HeuristicSearch::simplex, nb::call_guard<nb::gil_scoped_acquire>())
        .def("diagonal", &HeuristicSearch::diagonal, nb::call_guard<nb::gil_scoped_acquire>());
}