---
title: Building the Calculix adapter with PaStiX (and CalculiX 2.17)
permalink: adapter-calculix-pastix-build.html
keywords: adapter, calculix, pastix
summary: "Since version 2.17, CalculiX can be built with the PaStiX library to increase performance with CUDA. Building the preCICE adapter
is somehow harder with this version. This page gives a detailed walkthrough to build the modified adapter."
---

## Overview

To build the preCICE adapter for CalculiX with PaStiX, it is first necessary to build PaStiX. A manual build is necessary for most steps, as specific compilation flags must be used for the dependencies: pre-built packages cannot be used. In particular, PaStiX (and CalculiX) must be compiled with flags to use 8-byte integers. On this page, we provide a step-by-step building guide of the adapter, tested on Ubuntu 20.04.
This should work with some modifications on other systems. You may also need to tweak this depending on your needs, e.g. if you want to build some dependencies yourself instead of using a package manager.

{% note %}
This page was built following CalculiX' documentation with adaptations for preCICE and some additional tips from experience. Support of PaStiX is still experimental and feedback is highly valuable!
{% endnote %}

## Required packages

Make sure you have a working installation of preCICE. Also run these installation commands, (after a call to `sudo apt update` and `sudo apt upgrade`) :
`sudo apt install build-essential cmake git gfortran flex bison zlib1g-dev nvidia-cuda-toolkit-gcc libspooles-dev libyaml-cpp-dev`

{% note %}
You can get these elsewhere or build them from source. In particular, it is probably wise to have a more up to date CUDA installation than the one provided on the Ubuntu repository. Again, feedback is appreciated!
{% endnote %}

## Downloading CalculiX source

This guide assumes Calculix's source code is in the user's home folder `/home/user_name`, with the alias `~`. If you don't want to follow this convention, you may have to adapt slightly the instructions below. Download can be done on command line:

```bash
    cd ~ && wget http://www.dhondt.de/ccx_2.17.src.tar.bz2
    bunzip2 ccx_2.17.src.tar.bz2
    tar -xvf ccx_2.17.src.tar

```

This contains the source code of CalculiX, but also scripts useful for building some libraries below with the correct flags.

## Building PaStiX dependencies

It is assumed that all these libraries will be built on the user home folder, `~`. Minor modifications will be required otherwise.
PaStiX requires OpenBLAS, hwloc, Scotch and PaRSEC. All of these will be built before PaStiX itself. We need to build them from source because specific flags are necessary to work with CalculiX' version of PaStiX, mostly the use of 8 bytes integers.

### Building OpenBLAS

Clone OpenBLAS source code and build it with 8-byte integers.

```bash
    cd ~ 
    git clone https://github.com/xianyi/OpenBLAS.git 
    mv OpenBLAS ./OpenBLAS_i8
    cd OpenBLAS_i8 
    make -j 4 INTERFACE64=1 
    sudo make install

```

You can also specify a custom installation path (to avoid calling `sudo`) by using the `PREFIX=path/to/install` option of the Makefile.

### Building hwloc

This library will be put in a subfolder of the PaStiX folder.

```bash
    mkdir -p ~/PaStiX/ 
    cd ~/PaStiX/ 
    wget https://download.open-mpi.org/release/hwloc/v2.1/hwloc-2.1.0.tar.bz2
    bunzip2 hwloc-2.1.0.tar.bz2 && tar -xf hwloc-2.1.0.tar
    cp ~/CalculiX/ccx_2.17/src/make_hwloc.sh ~/PaStiX/hwloc-2.1.0/make_hwloc.sh
    cd hwloc-2.1.0
    ./configure --prefix=$HOME/PaStiX/hwloc_i8 CC=gcc CXX=g++
    make -j8
    make install

```

### Building PaRSEC

```bash
    cd ~/PaStiX && git clone -b pastix-6.0.2 --single-branch https://bitbucket.org/mfaverge/parsec.git
    cd parsec
    cp ~/CalculiX/ccx_2.17/src/make_parsec.sh ~/PaStiX/parsec
    ./make_parsec.sh

```

### Building Scotch

```bash
    cd ~/PaStiX
    wget https://gitlab.inria.fr/scotch/scotch/-/archive/master/scotch-master.tar.gz
    tar -xf scotch-master.tar.gz
    cp ~/CalculiX/ccx_2.17/src/make_scotch.sh ~/PaStiX/scotch-master
    cd scotch-master
    ./make_scotch.sh

```

## Building PaStiX

```bash
    cd ~/PaStiX
    git clone https://github.com/Dhondtguido/PaStiX4CalculiX 
    mv PaStiX4CalculiX pastix_src
    cd pastix_src
    rm make_pastix.sh
    cp ~/CalculiX/ccx_2.17/src/make_pastix.sh .
    ./make_pastix.sh

```

{% include warning.html content="The repository of PaStiX contains a `make_pastix.sh` file; be sure to use the one in the CalculiX folder and not the one from PaStiX!" %}

### Troubleshooting

- On some occasions, a Python script called by CMake (`cmake_modules/morse_cmake/modules/precision_generator/genDependencies.py`) generates an incorrect Makefile because of errors in regular expressions. This should be fixed by calling `pip install regex` (which requires installing the `python3-pip` Ubuntu package). You may also need to replace the `import re` line in that script by `import regex as re`, but the necessity seems to fluctuate among different machines.
- Some parts of the code require older versions of the GNU compilers. You may have to replace `gcc` by `gcc-7` and similarly for `g++` and `gfortran` in the `make_pastix.sh` script. This requires installing the relevant Ubuntu packages.

## Building ARPACK, a CalculiX dependency

Calculix relies on ARPACK, and when built with PaStiX, we cannot rely on standard distributions of that library, because it doesn't feature 8-byte integers by default. We need to build the [ARPACK code](https://en.wikipedia.org/wiki/ARPACK) (update: original page seems to be down) ourselves:

```bash
cd ~
wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz
wget https://web.archive.org/web/20220526222500fw_/https://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz
zcat arpack96.tar.gz | tar -xvf -
zcat patch.tar.gz    | tar -xvf -

```

Before building the library, the following modifications are required:

- In `ARmake.inc`, change `PLAT` by the appropriate suffix (the adapter's makefile assumes INTEL)
- In `ARmake.inc`, change the Fortran compilation flags to add `-fdefault-integer-8`. You may also need to remove the flag `-cg89`.
- If you extracted the archive on another folder than your home repository, update `home` in `ARmake.inc` accordingly.
- In the file `UTIL/second.f`, comment with a star the line containing `EXTERNAL ETIME`.

Once all of these are done, run `make lib` in the `ARPACK` folder.

## Building the adapter

To build the adapter, clone its repository and checkout the `2.17` branch :

```bash
    git clone -b v2.17 https://github.com/precice/calculix-adapter
    cd calculix-adapter
```

### Fixing the code

Due to some conflicts between CalculiX, PaStiX and the adapter (both CalculiX and PaStiX have a `pastix.h` file, and neither of them is local from the point of view of the adapter), some changes are required in the CalculiX codebase. We provide a script, `pastix_pre_build.sh` that does these changes. Run it.

```bash
    ./pastix_pre_build.sh
```

### Compilation

To build the adapter, use the provided `Makefile_i8_PaStiX`: the regular Makefile would build the adapter without PaStiX. Assuming you followed the previous steps, it should be useable without modifications other than giving Calculix' path; otherwise, some other paths updates could be required. You also need to ensure then Makefile finds the required dependencies when calling `pkg-config`. This can be done by changing the `PKG_PATH_CONFIG` environment variable. Assuming you used suggested paths, this would look like this:

```bash
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:~/PaStiX/pastix_i8/lib/pkgconfig/:~/PaStiX/hwloc_i8/lib/pkgconfig/:~/PaStiX/parsec_i8/lib/pkgconfig/
```

Then you can build the adapter with this command in the `calculix-adapter` (and be sure to checkout the `2.17` branch) folder :

```bash
    make -f Makefile_i8_PaStiX -j 4 CCX=~/CalculiX/ccx_2.17/src
```

Once the build is successful, the adapter should be in `./bin/ccx_preCICE`.

### Updating shared libraries

Running the adapter at this point should fail because of a missing shared library: `libparsec.so.2`. You can fix this by adding its path to the environment variable `LD_LIBRARY_PATH`: `export LD_LIBRARY_PATH=:$LD_LIBRARY_PATH:~/PaStiX/parsec_i8/lib`. You may have to do this everytime you load a new terminal, which is why we advise you to make this change permanent, for instance by adding this export commande at the end of your `.bashrc` file.
