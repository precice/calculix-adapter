# Automatic packaging of the adapter

This contains the necessary files to build a Debian package containing the adapter. (Built without PaStiX)

## Note for maintainers 

- Keep up to date the list of dependencies in the control file with `dpkg-shlibs`. Care must be taken when preCICE gets updated, typically.
- Update the changelog.
- The script `make_deb.sh` compresses the changelog and manpage and puts them at the right places. Building the package should be as simple as running this script.
    - Dependencies: `lintian` (to validate the package) and `pandoc` (to generate the `man` page).
