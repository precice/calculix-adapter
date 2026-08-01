---
title: CalculiX support
permalink: adapter-calculix-support.html
aliases:
  - /adapter-calculix-support.html
keywords: adapter, calculix, support, versions
summary: "Supported CalculiX versions and porting the CalculiX adapter to a different CalculiX version."
---

The CalculiX adapter directly modifies the source files of CalculiX and, as such, is made for a specific version of CalculiX.
This page includes some hints on porting the adapter to different versions.

{% tip %}
Have you upgraded the supported CalculiX version to a newer one? Contribute back to the repository and make it available for everyone!
{% endtip  %}

## Porting the adapter to a new CalculiX version

When a new CalculiX version is released, some files need to be updated:

- `CalculiX.h`
- `nonlingeo_precice.c`
- `cxx_2._version_.c` (version-specific: a new file replaces the previous one)

The adapter files in the `adapter/` folder should continue working independently of the CalculiX version.

Overall, the main loop of the solver must be modified to insert preCICE calls (time-stepping, communication, checkpointing, ...). The adapter files describe how these must be done (e.g. read data from CalculiX to write them into preCICE buffers, etc.).

{% tip %}
It is easier and safer to re-apply the (few) adapter-specific changes in the source files of the specific CalculiX version, instead of updating the source files and resolving the merge conflicts.
{% endtip %}

### File CalculiX.h

Simply copy the one from CalculiX and add the declaration of a function called `nonlingeo_precice`, which is the same as the existing `nonlingeo` but with two additional arguments at the end of the list: `char *preciceParticipantName, char *configFilename`.

### File nonlingeo_precice.c

Copy the file `nonlingeo.c` from CalculiX and add the adapter code. Use as reference the previous version: the difference between your `nonlingeo_precice.c` and the corresponding `nonlingeo.c` should be similar to the difference between the previous `nonlingeo_precice.c` and the `nonlingeo.c` of the previous CalculiX version. A good way to see what must be done is to use a graphical diff tool (e.g., [Meld](https://meldmerge.org/)). To make the diff smaller, try formatting the CalculiX code with the same `.clang-format` specification, or hiding whitespace changes and blank lines (e.g., with `diff -w -B`).

Adapter-specific changes are documented with a comment of the form `/* Adapter: doing XXX*/`. Finding these (i.e., looking for occurrences of the word `Adapter`) in the current version of the adapter can help ensure you didn't miss anything. *Return the courtesy by maintaining these comments!*

### File ccx_X.YY.c

Once again, copy the corresponding file, and add the adapter-specific changes, as above:

- Define the adapter-related variables
- Parse command line arguments to see if and how preCICE is being used
- If preCICE is being used, use the adapter-specific solver loop. This takes the form `if (preciceUsed) {custom loop} else if (some CCX code){...}`, replacing `if (some CCX code) {...}`.

## Bumping versions

There are several files where the CalculiX version or the adapter version and other details need to be updated. See the pull request template `.github/PULL_REQUEST_TEMPLATE/release.md` for a checklist.
