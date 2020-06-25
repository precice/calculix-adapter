cmake_minimum_required(VERSION 3.10)

include(ExternalProject)

project(CXXAdapter CXX C Fortran)

find_package(Threads REQUIRED)
find_package(precice REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

set(CCX_VERSION "2.16")

function(download_once url file)
  if(NOT EXISTS ${file})
    file(DOWNLOAD ${url} ${file} SHOW_PROGRESS)
  endif()
endfunction(download_once)

#
# Setup YAMLcpp
#

message(STATUS "Downloading YAML-cpp")
download_once(https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.2.zip external/yamlcpp/yaml-cpp.zip)
externalproject_add(
  YAMLCpp
  PREFIX external/yamlcpp
  URL ${CMAKE_BINARY_DIR}/external/yamlcpp/yaml-cpp.zip
  CMAKE_ARGS -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF -DYAML_CPP_BUILD_TOOLS=OFF -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND ""
)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/yamlcpp/src/YAMLCpp/include)

add_library(yamlcpp STATIC IMPORTED)
set_target_properties(
  yamlcpp PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
  ${CMAKE_BINARY_DIR}/external/yamlcpp/src/YAMLCpp/include
  IMPORTED_LOCATION
  ${CMAKE_BINARY_DIR}/external/yamlcpp/src/YAMLCpp/libyaml-cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
)

#
# Setup Spooles
#

message(STATUS "Downloading Spooles")
download_once(http://www.netlib.org/linalg/spooles/spooles.2.2.tgz external/spooles/spooles.tgz)
configure_file(patch/Spooles/Make.inc.in external/spooles/patch/Make.inc @ONLY)
externalproject_add(
  Spooles
  PREFIX external/spooles
  URL ${CMAKE_BINARY_DIR}/external/spooles/spooles.tgz
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/external/spooles/patch/Make.inc ${CMAKE_BINARY_DIR}/external/spooles/src/Spooles/Make.inc
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make lib
  BUILD_IN_SOURCE True
  INSTALL_COMMAND ""
)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/spooles/src/Spooles)

add_library(spooles STATIC IMPORTED)
set_target_properties(
  spooles PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
  ${CMAKE_BINARY_DIR}/external/spooles/src/Spooles
  IMPORTED_LOCATION
  ${CMAKE_BINARY_DIR}/external/spooles/src/Spooles/spooles${CMAKE_STATIC_LIBRARY_SUFFIX}
)

#
# Setup ARPACK
#

message(STATUS "Downloading ARPACK")
download_once(https://github.com/opencollab/arpack-ng/archive/master.zip external/arpack/arpack.zip)
externalproject_add(
  ARPACK
  PREFIX external/arpack
  URL ${CMAKE_BINARY_DIR}/external/arpack/arpack.zip
  CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release
  BUILD_COMMAND ${CMAKE_COMMAND} --build . --target arpack
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND ""
)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/arpack/src/ARPACK)

add_library(arpack STATIC IMPORTED)
set_target_properties(
  arpack PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES
  ${CMAKE_BINARY_DIR}/external/arpack/src/ARPACK
  IMPORTED_LOCATION
  ${CMAKE_BINARY_DIR}/external/arpack/src/ARPACK/libarpack${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_LINK_LIBRARIES
  "${LAPACK_LIBRARIES}"
)

#
# Setup CCX
#

set(CCX_ROOT ${CMAKE_BINARY_DIR}/ccx)
set(CCX_SOURCE_FILE ${CCX_ROOT}/ccx_source.tar.bz2)
set(CCX_SOURCE_DIR ${CCX_ROOT}/CalculiX/ccx_${CCX_VERSION}/src)
message(STATUS "Downloading CalculiX ${CCX_VERSION}")
download_once(http://www.dhondt.de/ccx_${CCX_VERSION}.src.tar.bz2 ${CCX_SOURCE_FILE})
if(NOT IS_DIRECTORY ${CCX_ROOT}/CalculiX)
  message(STATUS "Preparing CalculiX sources")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xz ${CCX_SOURCE_FILE}
    WORKING_DIRECTORY ${CCX_ROOT}
  )
  file(RENAME ${CCX_SOURCE_DIR}/ccx_${CCX_VERSION}.c ${CCX_SOURCE_DIR}/ccx_${CCX_VERSION}.c.orig)
  # file(RENAME ${CCX_SOURCE_DIR}/CalculiX.h ${CCX_SOURCE_DIR}/CalculiX.h.orig)
endif()

message(STATUS "Detecting CalculiX files")
# Some source files only work if they are included by others.
# So we have to use the sources defined in the Makefile.inc
file(READ "${CCX_SOURCE_DIR}/Makefile.inc" _ccx_sources)
string(REGEX MATCHALL "[A-Za-z0-9_-]+.(f|cpp|c)" _ccx_local_sources "${_ccx_sources}")
set(CCX_SOURCES "")
foreach(file ${_ccx_local_sources})
  # CXX files are not required for the build
  if(NOT "${file}" MATCHES ".cpp$")
    set(CCX_SOURCES ${CCX_SOURCES} ${CCX_SOURCE_DIR}/${file})
  endif()
endforeach(file)
message(STATUS "Detecting adapter files")
file(GLOB ADAPTER_SOURCES "*.c" "*.f" "adapter/*")

message(STATUS "Setting up preCICE adapter")
add_executable(ccx_precice ${CCX_SOURCES} ${ADAPTER_SOURCES})
add_dependencies(ccx_precice ARPACK YAMLCpp Spooles)
set_target_properties(
  ccx_precice
  PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED TRUE
  CXX_EXTENSIONS FALSE
)
target_include_directories(
  ccx_precice PRIVATE
  ${CMAKE_SOURCE_DIR}
  ${CCX_SOURCE_DIR}
)
target_link_libraries(
  ccx_precice PRIVATE
  OpenMP::OpenMP_C
  Threads::Threads
  precice::precice
  spooles
  arpack
  yamlcpp
)
target_compile_options(
  ccx_precice PRIVATE
  $<$<COMPILE_LANGUAGE:C>:-fcommon>
  $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
)
target_compile_definitions(
  ccx_precice PRIVATE
  ARCH=Linux
  SPOOLES
  ARPACK
  MATRIXSTORAGE
)
