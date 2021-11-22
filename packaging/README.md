# Automatic packaging of the adapter

This contains the necessary files to build a Debian package containing the adapter. (Built without PaStiX)

## Note for maintainers 

- Keep up to date the list of dependencies in the control file with `dpkg-shlibs`. In theory, it will detect preCICE with the version used for the build, so we update it manually to the lowest possible version. For instance, if the adapter is built with preCICE 2.3.0, we update the dependency as >= 2.0.0
- Update the changelog.
- The script `make_deb.sh` compresses the changelog and manpage and puts them at the right places. Building the package should be as simple as running this script.
    - Dependencies: `lintian` (to validate the package) and `pandoc` (to generate the `man` page).
- To update CalculiX version: update the Github action workflow file (in `.github/workflows`). Update the `PACKAGE_NAME` variable, the paths (e.g. replace `wget http://www.dhondt.de/ccx_2.17.src.tar.bz2` by `wget http://www.dhondt.de/ccx_2.XX.src.tar.bz2` with XX the corresponding minor version of CalculiX and folder paths.
