# Updating the CalculiX adapter: some guidelines

## Versioning scheme

The adapter has a versioning scheme inherited from CalculiX. It is of the form `CCX_MAJOR.CCX_MINOR.ADAPTER_PATCH`. For instance, the first release was 2.19.0 because it was based on CalculiX 2.19. If the adapter gets a new release before a new CalculiX version arrives, it shall be named with a bump in the patch number (e.g. 2.19.1), or with the new CalculiX version and patch number 0 otherwise. This would be 2.20.0, in this example. New features should be merged on the `develop` branch and get merged in `master` at release time.

## Checklist for releasing an adapter update

Whether or not there is a CalculiX update, the following files need an update (typically update the value of a global variable)

- The Github workflow file for building Debian Packages: `.github/workflows/ubuntu_build.yml`
- Files related to the build of the Debian Package, in the `packaging` folder: `calculix_package/DEBIAN/control`, `make_deb.sh` and `changelog.Debian` (where an entry should be added)
- On the preCICE website, the variable `calculix_adapter_version` in the `_config.yml` file
- The repository's `README`.

## Checklist for upgrading to a newer CalculiX version

Apart from updating the CalculiX file in the source code (see below) and the things mentioned above, the following should be updated:

- The `CCX_VERSION` in the Makefile(s)
- The `pastix_pre_build.sh` file
- On the preCICE website, the variable `calculix_version` in the `_config.yml` file

### How to port the adapter to a new CalculiX version

When a new CalculiX version is released (at least for 2.xx), some files need to be updated:

- The header `CalculiX.h`
- `nonlingeo_precice.c`
- A new file `cxx_2._version_.c` must replace the old one

On the other hand, adapter files in the `adapter/` folder should be working independantly of the current CalculiX version. Roughly speaking, the files to change contain the main loop of the solver and must be modified to insert preCICE calls (time-stepping, communication, checkpointing, ...) whereas the adapter files describe how these must be done (e.g. read data from CalculiX to write them into preCICE buffers, etc.).
To make sure we get the correct CalculiX code, these 3 files must be rewritten for each CalculiX version. For instance, the file `ccx_2.19.c` is not an update of `ccx_2.17.c`, but a copy of the `ccx_2.19.c` of CalculiX source code, modified with adapter code inserted.

#### `CalculiX.h`

Simply copy the one from CalculiX and add the declaration of a function called `nonlingeo_precice` which is the same as the existing `nonlingeo` but with two additional arguments: `char *preciceParticipantName, char *configFilename`. It is better to add them at the end.

#### `nonlingeo_precice.c`

Copy the file `nonlingeo.c` from CalculiX and add adapter code. Use as reference the previous version: the difference between your `nonlingeo_precice.c` and the corresponding `nonlingeo.c` should be similar to the difference between the previous `nonlingeo_precice.c` and the `nonlingeo.c` of the previous CalculiX version. A good way to see what must be done is to use the `diff -w -B` command to get these differences. (You may have to format the CalculiX source code to reduce noise.)
As of now, all changes are documented with a comment of the form `/* Adapter: doing XXX*/`. Finding these (i.e. looking for all occurences of the word `Adapter`) in previous version of the adapter can help ensure you didn't miss anything. **Make sure to return the courtesy by keeping these comments.**

#### 'ccx_X.YY.c`

Once again, copy the corresponding file and add preCICE code, as above. As of now, this consists of 3 things:

- Define some preCICE-related variables like the participant name.
- Parse command line arguments to see if and how preCICE is used
- If preCICE is used, use our custom solver loop instead of ones from CalculiX. This takes the form `if (preciceUsed) {custom loop} else if (some CCX code){...}` replacing `if (some CCX code) {...}`.