---
title: Troubleshooting the CalculiX adapter
permalink: adapter-calculix-troubleshooting.html
aliases:
  - /adapter-calculix-troubleshooting.html
keywords: adapter, calculix, error
summary: "While working with the CalculiX adapter, you may run onto common issues. This is a collection of what we know could go wrong."
---

This list is definitely not complete. If after reading this, you still have issues, please [ask in the preCICE forum](https://precice.discourse.group/).

## Things to check

* Are you using the same version of CalculiX and of the CalculiX adapter? The adapter installation works by replacing files of the original CalculiX, so they should be compatible.
* Can you manually run the `ccx_preCICE` binary?
  * It should be in your `$PATH`
  * If autocompletion does not work (e.g. `ccx_` + TAB key), then it is probably not in your `$PATH`.
* Our tutorials also require CGX (pre- and post-processor of CalculiX).
  * Is CGX installed?
  * Is OpenGL (required by CGX) installed?

## Compiling with GCC 10 or newer

If you compile the adapter (v2.20.1 or earlier) with GCC 10 or newer, you will get the following error, originating from CalculiX:

```text
Error: Rank mismatch between actual argument at (1) and actual argument at (2) (rank-1 and scalar)
```

To work around this, you need to add `-fallow-argument-mismatch` to the `FFLAGS` inside `Makefile`:

```diff
- FFLAGS = -Wall -O3 -fopenmp $(INCLUDES)
+ FFLAGS = -Wall -O3 -fopenmp -fallow-argument-mismatch $(INCLUDES)
```

## Compiling with Intel OneAPI

If you compile the adapter with Intel OneAPI, you will get the following error at link time, originating from CalculiX:

```text
undefined reference to MAIN__
```

To work around this, you need to add the `-nofor-main` to the `FFLAGS` inside `Makefile`. The easiest way is to define `ADDITIONAL_FFLAGS`:

```shell
ADDITIONAL_FFLAGS="-nofor-main` make
```

to specify that the main program is not written in Fortran. Read more in the [Intel compiler reference](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2026-1/nofor-main.html).

## Undefined references to SPOOLES functions

On systems where both `spoolesMT.a` (multithreading) and `spooles.a` are available, you might get undefined references to SPOOLES. In that case, you might need to modify `SPOOLES_LIBS` in the `Makefile` to include the paths to both `spoolesMT.a` and to `spooles.a`, in this order. See the [`Makefile_MT` file of CalculiX](https://github.com/g0mb4/CalculiX/blob/58c684679097132eda8bb827f186b5feeaa1281f/src/Makefile_MT).
