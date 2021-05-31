---
title: Get CalculiX
permalink: adapter-calculix-get-calculix.html
keywords: adapter, calculix, building, spooles, arpack, yaml-cpp
summary: "Building CalculiX itself can already be quite a challenges. That's why we collected here some recipe."
---

The CalculiX adapter for preCICE directly modifies the source code of CalculiX and produces an alternative executable `ccx_preCICE`. Therefore, we first need to get and (optionally) build CalculiX from source.

[CalculiX](http://www.dhondt.de) consists of the solver, called "CCX" and a pre- and postprocessing software with graphical user interface "CGX".

- The installation procedure of CCX is described in its `src/README.INSTALL` files, but we also give a summary here.
- We don't modify CGX, so you can simply get a binary package (if needed, e.g. as a preprocessor in our FSI tutorials)

You don't need to build the "vanilla" CalculiX before building the adapter. But you do need to get all the dependencies and the source code of CCX.

## Dependencies

CalculiX itself depends on [SPOOLES2.2](http://www.netlib.org/linalg/spooles/spooles.2.2.html) and [ARPACK](https://www.caam.rice.edu/software/ARPACK/).

Additionally, our adapter also depends on [yaml-cpp](https://github.com/jbeder/yaml-cpp).

These can be found in many distributions as binary packages. For example, in Ubuntu, do:

```bash
sudo apt install libarpack2-dev libspooles-dev libyaml-cpp-dev
```

For example, in Arch or Manjaro, install `arpack` and `yaml-cpp`, and compile `spooles` using an AUR helper (e.g. `yay`):

```bash
sudo pacman -S arpack yaml-cpp
yay spooles
```

If `spooles` compilation breaks with `-Werror=format-security`, replace the flag with `-Wformat-security` in `CFLAGS` (file `/etc/makepkg.conf`).

### Building Spooles from source

<details markdown="1"><summary>If you cannot get a binary for Spooles, try these instructions.</summary>

Download SPOOLES, e.g:

```bash
wget http://www.netlib.org/linalg/spooles/spooles.2.2.tgz 
```

Extract it in a separate directory

```bash
mkdir SPOOLES.2.2
tar zxvf spooles.2.2.tgz -C SPOOLES.2.2
cd SPOOLES.2.2
```

Edit by hand configuration file `Make.inc` to change the compiler version in line 14-15

```make
CC = gcc
#CC = /usr/lang-4.0/bin/cc 
```

Now build the library:

```bash
make lib 
```

</details>

### Building ARPACK from source

<details markdown="1"><summary>If you cannot get a binary for ARPACK, try these instructions.</summary>

Download Arpack and patch:

```bash
wget https://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz 
wget https://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz 
```

Unpack them (they will be unpacked in the newly created directory `ARPACK`)

```bash
tar xzfv arpack96.tar.gz 
tar xzfv patch.tar.gz 
cd ARPACK
```

Edit by hand `ARmake.inc` to specify build instructions. The following changes will depend on the directory structure of your system:

- **Line 28**:  Change `home = $(HOME)/ARPACK` to directory where ARPACK in extracted
- **Line 115**: Change `MAKE    = /bin/make` to e.g. `MAKE  =  make` (if needed)
- **Line 120**: Change `SHELL   = /bin/sh` to e.g. `SHELL =  sh`  (if needed)
- **Lines 104 - 105**: Specify your fortran compiler and compiler flags, e.g. for the gnu systems:

```make
FC = gfortran
#FFLAGS = -O -cg89 
```

- **Line 35**: Modify the platform suffix for the library and remember it, since Calculix adapter makefile will depend on it ( by default it will use suffix INTEL for Linux and MAC for mac systems). For example change
`PLAT = SUN4` to `PLAT = INTEL`
- You will probably get linking errors related to ETIME, which you can bypass: In the  file `UTIL/second.f` append `*` to the beginning of line 24 ( that comments it out )

    ```fortran
    *        EXTERNAL           ETIME
    ```

Now we are ready to build the library with `make lib`
</details>

### Building yaml-cpp from source

<details markdown="1"><summary>If you cannot get a binary for yaml-cpp, try these instructions.</summary>

Get the latest version of yaml-cpp and build it as a shared library. For example:

```bash
wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.2.zip
unzip yaml-cpp-0.6.2.zip
cd yaml-cpp-yaml-cpp-0.6.2
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=ON ..
make
```

After building, make sure that you make yaml-cpp discoverable by setting e.g. your `LD_LIBRARY_PATH`. You don't need this for the CalculiX adapter, but you would need it e.g. for the OpenFOAM adapter.

**Note**: If you use Boost 1.67 or newer, then you also need to install yaml-cpp 0.6 or newer. Similarly, for an older Boost version, you also need an older yaml-cpp. Unfortunately, this is not related to the adapter's code.
</details>

## Building CalculiX with the preCICE adapter

### Get the source

Once the libraries are installed, you can finally install Calculix with preCICE adapter (adapt the `VERSION` in the link, see the beginning of the [adapter's README.md](https://github.com/precice/calculix-adapter/blob/master/README.md) to find out which one you need).

```bash
cd ~
wget http://www.dhondt.de/ccx_VERSION.src.tar.bz2
tar xvjf ccx_VERSION.src.tar.bz2 
```

The source code is now in the `~/CalculiX/ccx_VERSION/src` directory. The adapter's [`Makefile`](https://github.com/precice/calculix-adapter/blob/master/Makefile) is looking for CCX in this directory by default, so modify it if needed.

### Building the "vanilla" CalculiX (optional)

If you want to build the "vanilla" (i.e. without preCICE) CalculiX, you can now run `make` inside the `src/` directory. Depending on how you installed the dependencies above (using `apt` or from source), you might get compilation errors, such as `spooles.h:26:10: fatal error: misc.h: No such file or directory`. Often these errors can be easily fixed by modifying CalculiX `Makefile`. Please refer to [our adapter's makefile options](adapter-calculix-get-adapter.html#makefile-options) for a list of library and include flag you might have to set depending on your installation procedure.

### Building the modified CalculiX

Continue to the page [get the adapter](adapter-calculix-get-adapter.html).
