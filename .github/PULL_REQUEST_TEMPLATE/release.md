## Release checklist

I updated the adapter version in the following:

- [ ] The version listed in `README.md`
- [ ] The `CCX_VERSION` in the Makefile(s)
- [ ] The `pastix_pre_build.sh` file
- [ ] Files related to the build of the Debian Package, in the `packaging` folder:
  - `debian/DEBIAN/control`
  - `make_deb.sh`
  - `changelog.Debian` (where an entry should be added)
- [ ] The Github workflow file:
  - `.github/workflows/ubuntu_build.yml`
  - `.github/workflows/request-deb-package.yaml`
- [ ] The documentation pages
- [ ] `CHANGELOG.md`

Outside this repository:

- [ ] [tutorials](https://github.com/precice/tutorials/tree/develop/tools/tests): Update the default CalculiX version in the system tests:
  - `components.yaml`
  - `reference_versions.yaml` (if the reference results need to be updated)

System tests:

- [ ] I triggered the system tests by adding the `trigger-system-tests` label.
