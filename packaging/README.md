# Automatic packaging of the adapter

This contains the necessary files to build a Debian package containing the adapter. (Built without PaStiX)

## Note for maintainers 

- Keep up to date the list of dependencies in the control file with `dpkg-shlibs`. Care must be taken when preCICE gets updated, typically.
- Update the changelog.
- The script `make_deb.sh` does the required actions to compress the changelog and manpage and put them at the right place. Building the package should be as simple as running it.
- `lintian` and `pandoc` are required to run that script. The first checks the validity of the package and the second compiles the `man` page from Markdown into the appropraite format
