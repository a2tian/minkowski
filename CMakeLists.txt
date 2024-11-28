# Set the minimum CMake version and policies for highest tested version
cmake_minimum_required(VERSION 3.15...3.27)

# Set up the project and ensure there is a working C++ compiler
project(minkowski LANGUAGES CXX)

# Warn if the user invokes CMake directly
if (NOT SKBUILD)
message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
  
find_package(nanobind CONFIG REQUIRED)

find_library(NTL_LIB ntl)
find_library(GMP_LIB gmp)
find_library(M_LIB m)
find_package(OpenMP)
find_package(Threads)

file(GLOB SOURCES "cc/*.cc")

nanobind_build_library(nanobind-static)
nanobind_add_module(_minkowski NB_STATIC ${SOURCES})
target_link_libraries(_minkowski PRIVATE Threads::Threads OpenMP::OpenMP_CXX ${NTL_LIB} ${GMP_LIB} ${M_LIB})

nanobind_add_stub(_minkowski_stub 
  MODULE _minkowski
  OUTPUT _minkowski.pyi
  MARKER_FILE py.typed
  PYTHON_PATH $<TARGET_FILE_DIR:_minkowski>
  DEPENDS _minkowski
)  

install(TARGETS _minkowski LIBRARY DESTINATION minkowski)
