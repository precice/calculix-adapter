---
title: Get the CalculiX adapter
permalink: adapter-calculix-get-adapter.html
keywords: adapter, calculix, building
summary: "The CalculiX adapter provides the executable `ccx_preCICE`. You can get the adapter either from a Debian package (on Ubuntu), or build it from source."
---

After [installing preCICE](installation-overview.html) and [getting the CalculiX source and the required dependencies](adapter-calculix-get-calculix.html), you can now build the adapter, i.e. a modified CCX executable.

There are two ways to get the adapter: (a) get a binary package (Ubuntu-only), or (b) build it from source. The latest adapter version is v{{ site.calculix_adapter_version }} and the versioning format is `<CalculiX major.minor version>.<adapter revision>`.

## Get a binary package

You can download version-specific Ubuntu (Debian) packages from each [adapter release](https://github.com/precice/calculix-adapter/releases/latest). To install, open it in your software center.

Alternatively, download & install it from the command line. For Ubuntu 22.04 (jammy):

```bash
wget https://github.com/precice/calculix-adapter/releases/download/v{{ site.calculix_adapter_version }}/calculix-precice2_{{ site.calculix_adapter_version }}-1_amd64_jammy.deb
sudo apt install ./calculix-precice2_{{ site.calculix_adapter_version }}-1_amd64_jammy.deb
```

This requires that also preCICE itself has been installed from a Debian package.

{% tip %}
Replace `jammy` with `focal` to get the package for Ubuntu 20.04.
{% endtip  %}

{% note %}
We started offering Debian packages for the CalculiX adapter since v2.19.0. Please [give us your feedback](community-channels.html)!
{% endnote %}

## Building the adapted CalculiX

1. Download and unzip the latest state of the adapter (e.g. in the `CalculiX` folder), currently supporting CalculiX v{{site.calculix_version}}:

    ```bash
    wget https://github.com/precice/calculix-adapter/archive/refs/heads/master.tar.gz
    tar -xzf master.tar.gz
    cd calculix-adapter-master
    ```

2. Edit the `Makefile` to set the paths to dependencies.
   - If you have the CalculiX source in `~/CalculiX/` and the dependencies in your global paths, you don't need to change anything.
   - Otherwise, set `CCX` and, if built from source, the include and lib flags for SPOOLES, ARPACK, and yaml-cpp.
   - Read below if you are [using GCC 10 or later](#compiling-with-gcc-10-or-newer).
3. Clean any previous builds with `make clean`.
4. Build with `make` (e.g. `make -j 4` for parallel).
5. You should now have a new executable `ccx_preCICE` in the `bin/` folder of the adapter. You may move this file to a path known by your system, or [add this to your `PATH`](https://unix.stackexchange.com/a/26059/36693) (careful when doing this!).

### Building the adapter with PaStiX

Since version 2.17 of CalculiX, it is possible to link the PaStiX solver for increased performance, mostly through GPUs. Building the adapter with PaStiX is quite tedious, as most dependencies of PaStiX and PaStiX itself must be built from source. Check our [detailed instructions on building the adapter with PaStiX](adapter-calculix-pastix-build.html).

### Makefile options

The adapter is built using GNU Make. The `Makefile` contains a few variables on top, which need to be adapted to your system:

 1. `CCX`: Location of the original CalculiX solver (CCX) source code ("src" directory)
    - Example: `$(HOME)/CalculiX/ccx_{{site.calculix_version}}/src`
 2. `SPOOLES_INCLUDE`: Include flags for SPOOLES
    - Example 1: `SPOOLES_INCLUDE   = -I/usr/include/spooles/` (installed)
    - Example 2: `SPOOLES_INCLUDE   = -I$(HOME)/SPOOLES.2.2/` (source)
 3. `SPOOLES_LIBS`: Library flags for SPOOLES
    - Example 1: `SPOOLES_LIBS      = -lspooles` (installed)
    - Example 2: `SPOOLES_LIBS      = $(HOME)/SPOOLES.2.2/spooles.a` (source)
 4. `ARPACK_INCLUDE`: Include flags for ARPACK
    - Example 1: `ARPACK_INCLUDE    =` (installed, nothing needed)
    - Example 2: `ARPACK_INCLUDE    = -I$(HOME)/ARPACK` (source)
 5. `ARPACK_LIBS`: Library flags for ARPACK
    - Example 1: `ARPACK_LIBS       = -larpack -llapack -lblas` (installed)
    - Example 2: `ARPACK_LIBS       = $(HOME)/ARPACK/libarpack_INTEL.a` (source)
 6. `YAML_INCLUDE`: Include flags for yaml-cpp
    - Example 1: `YAML_INCLUDE      = -I/usr/include/` (installed)
    - Example 2: `YAML_INCLUDE      = -I$(HOME)/yaml-cpp/include` (source)
 7. `YAML_LIBS`: Library flags for yaml-cpp
    - Example 1: `YAML_LIBS         = -lyaml-cpp` (installed)
    - Example 2: `YAML_LIBS         = -L$(HOME)/yaml-cpp/build -lyaml-cpp` (source)

You may also want to adjust the compiler `FC` from `mpifort` to `mpif90` or to any other compiler your system uses.

### Compiling with GCC 10 or newer

If you compile with GCC 10 or newer, you will get the following error, originating from CalculiX:

```text
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
```

To work around this, you need to add `-fallow-argument-mismatch` to the `FFLAGS` inside `Makefile`:

```diff
- FFLAGS = -Wall -O3 -fopenmp $(INCLUDES)
+ FFLAGS = -Wall -O3 -fopenmp -fallow-argument-mismatch $(INCLUDES)
```

### Notes on preCICE versions

<details markdown="1"><summary>In case you are using some very old preCICE version, please upgrade. Our <a href="https://precice.discourse.group/" title="preCICE forum">community</a> is happy to help you. Click here and keep reading if you loved preCICE v1.x and (optionally) wish The Beatles were still around.</summary>

1. This adapter expects the preCICE C bindings in `[prefix]/include/precice/SolverInterfaceC.h` and gets this path from pkg-config. In other words, this assumes that preCICE (at least v1.4.0) has been built & installed with CMake (e.g. using a Debian package). In case you want to keep using preCICE built with SCons, see the changes invoked by [Pull Request #14](https://github.com/precice/calculix-adapter/pull/14).
2. Starting from preCICE v1.2.0, the name (and the respective paths) of the language "adapters" have changed to language "bindings". This affects the line `#include "precice/bindings/c/SolverInterfaceC.h"` in `calculix-adapter/adapter/PreciceInterface.c`. To compile with older preCICE versions, change `bindings` to `adapters`.

</details>
