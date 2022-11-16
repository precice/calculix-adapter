---
title: Building the CalculiX adapter on SuperMUC
permalink: adapter-calculix-supermuc.html
keywords: adapter, calculix, cluster, modules
summary: "This page explains how to build the CalculiX adapter on SuperMUC. Even though SuperMUC was shut down in 2019, this page may still be useful for other clusters."
---

{% warning %}
This page needs updates for preCICE v2.
{% endwarning %}

In order to install CalculiX and the adapter on superMUC, a number of dependencies are first required. Initially, [preCICE must be installed](installation-overview.html)

Additionally, [SPOOLES](http://www.netlib.org/linalg/spooles/spooles.2.2.html), [ARPACK](https://en.wikipedia.org/wiki/ARPACK) and [yaml-cpp](https://github.com/jbeder/yaml-cpp) are required.

To install SPOOLES, some changes are necessary.

1. `makefile`: `~/SPOOLES.2.2/Tree/src/makeGlobalLib` contains an error: file `drawTree.c` does not exist and should be replaced by `tree.c`.
2. Changes to the `Make.inc` file must be done according to the CalculiX install [Manual](http://www.dhondt.de/INST_CCX_2_8_MAC_02_10_2015.pdf), page 16 and 17.

In installing ARPACK, the HOME directory needs to be specified in the "ARmake.inc" file. No changes are necessary for the Makefile. To install ARPACK, run "make lib" in the ARPACk directory.

Any problems with the installation of SPOOLES and ARPACK can be searched in the installation [instructions](http://www.dhondt.de/INST_CCX_2_8_MAC_02_10_2015.pdf).

To install yaml-cpp, run in the source directory:

```bash
mkdir build
cd build
cmake ..

make
make install
```

yaml-cpp 0.5.3 is known to work. Newer version may also work. yaml-cpp can be downloaded from

```bash
wget https://github.com/jbeder/yaml-cpp/archive/release-0.5.3.tar.gz -O - | tar xz 
```

## Module List

The following modules available in superMUC are known to work for the CalculiX adapter installation.

1. python/3.5_anaconda_nompi
2. scons/3.0.1  
3. valgrind/3.10
4. petsc/3.8
5. boost/1.65_gcc
6. gcc/6
7. mpi.intel/2017

## Makefile Changes

The paths to the CalculiX CCX, SPOOLES, ARPACK and YAML must be specified.
Line 61: "FC = mpifort" can be commented out and replaced with "FC = gfortran".

The path to the pkgconfig file needs to be stated. The command "export PKG_CONFIG_PATH=/path/to/lib/pkgconfig" must be provided. It is easier to install preCICE with the "CMAKE_INSTALL_PREFIX" set to the path where preCICE is installed.

## Adapter Installation

To install the adapter, the command with the following configurations is known to work:

```bash
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
```
