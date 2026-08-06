---
title: Get the CalculiX adapter
permalink: adapter-calculix-get-adapter.html
aliases:
  - /adapter-calculix-get-adapter.html
keywords: adapter, calculix, building
summary: "The CalculiX adapter provides the executable `ccx_preCICE`. You can get the adapter either from a Debian package (on Ubuntu), or build it from source."
---

After [installing preCICE](https://precice.org/installation-overview.html) and [getting the CalculiX source and the required dependencies](adapter-calculix-get-calculix.html), you can now build the adapter, i.e., a modified CalculiX executable.

There are two ways to get the adapter: (a) get a binary package (Ubuntu-only) or (b) build it from source. The adapter follows the versioning format `<CalculiX major.minor version>.<adapter revision>`.

## Get a binary package

You can download version-specific Ubuntu (Debian) packages from each [adapter release](https://github.com/precice/calculix-adapter/releases/latest). To install, open it in your software center.

Alternatively, download & install it from the command line. For Ubuntu 26.04 (Resolute Raccoon):

```bash
wget https://github.com/precice/calculix-adapter/releases/download/v2.20.2/calculix-precice3_2.20.2-1_amd64_resolute.deb
sudo apt install ./calculix-precice3_2.20.2-1_amd64_resolute.deb
```

This requires that also preCICE itself has been installed from a Debian package.

{% important %}
Replace `resolute` with `noble` to get the package for Ubuntu 24.04 (Noble Numbat), or with `jammy` for Ubuntu 22.04 (Jammy Jellyfish).
{% endimportant  %}

## Building the adapted CalculiX

1. Download and unzip the latest state of the adapter (e.g. in the `CalculiX` folder), and see the `README.md` for the supported CalculiX version:

    ```bash
    wget https://github.com/precice/calculix-adapter/archive/refs/heads/master.tar.gz
    tar -xzf master.tar.gz
    cd calculix-adapter-master
    ```

2. Edit the `Makefile` to set the paths to dependencies.
   - If you have the CalculiX source in `~/CalculiX/` and the dependencies in your global paths, you don't need to change anything.
   - Otherwise, set `CCX` and the include and lib flags for the dependencies.
3. Clean any previous builds with `make clean`.
4. Build with `make` (e.g., `make -j 4` for parallel).
5. You should now have a new executable `ccx_preCICE` in the `bin/` folder of the adapter. You may move this file to a path known by your system, or [add this to your `PATH`](https://unix.stackexchange.com/a/26059/36693) (careful when doing this!).

### Building the adapter with PaStiX

CalculiX can link to the PaStiX solver for increased performance using GPUs. Building the adapter with PaStiX is quite tedious, as most dependencies of PaStiX and PaStiX itself must be built from source. Check some [instructions on building the adapter with PaStiX](adapter-calculix-pastix-build.html).

### Makefile options

The adapter is built using GNU Make. The `Makefile` contains a few variables on top, which need to be adapted to your system:

 1. `CCX`: Location of the original CalculiX solver (CCX) source code ("src" directory)
    - Example: `$(HOME)/CalculiX/ccx_2.20/src`
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

{% version %}
The variables `YAML_INCLUDE` and `YAML_LIBS` are only relevant up to the adapter v2.20.1.
{% endversion %}

You may also want to adjust the compiler `FC` from `mpifort` to `mpif90` or to any other compiler your system uses.

See also the [troubleshooting](adapter-calculix-troubleshooting.html) page for further known issues.
