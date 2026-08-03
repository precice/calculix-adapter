# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

<!-- markdownlint-configure-file {"MD024": { "siblings_only": true } } -->

## [Unreleased]

### Added

- Added preliminary functionality to support volumetric coupling by extracting integration points from an element mesh ([#146](https://github.com/precice/calculix-adapter/pull/146)).
- Extended velocities writing for the 2D3D case ([#137](https://github.com/precice/calculix-adapter/pull/137)).
- Added Ubuntu 24.04 to the CI ([#140](https://github.com/precice/calculix-adapter/pull/140)).

### Changed

- Moved `struct SimulationData` in `nonlingeo_precice.h` to just before `Precice_Setup` ([#154](https://github.com/precice/calculix-adapter/pull/154))
- Improved Makefile with respect to compiler flags and version-aware flags ([#139](https://github.com/precice/calculix-adapter/pull/139)).
- Migrated to header-only yaml-cpp ([#143](https://github.com/precice/calculix-adapter/pull/143)).

### Fixed

### Removed

- Removed Ubuntu 24.04 from the CI ([#140](https://github.com/precice/calculix-adapter/pull/140)).

## [v2.20.1] 2024-03-20

### Added

- Added pre-commit hook config file for clang-format ([#117](https://github.com/precice/calculix-adapter/pull/117)).
- Added includes in more build steps ([#119](https://github.com/precice/calculix-adapter/pull/119)).

### Changed

- Ported the CalculiX adapter to preCICE v3 ([#114](https://github.com/precice/calculix-adapter/pull/114)).

## [v2.20.0] 2022-11-17

### Added

- Added modal dynamic simulations ([#91](https://github.com/precice/calculix-adapter/pull/91)), including support for implicit coupling ([#99](https://github.com/precice/calculix-adapter/pull/99)) and subcycling ([#105](https://github.com/precice/calculix-adapter/pull/105)).
- Added support for hexaedral elements in face meshes. Either hexaedral or tetrahedral elements can be used, but not a combination ([#91](https://github.com/precice/calculix-adapter/pull/91)).
- Added support for reading pressure ([#91](https://github.com/precice/calculix-adapter/pull/91)).
- Added the possibility of using a static step before the actual coupled simulation ([#101](https://github.com/precice/calculix-adapter/pull/101)).

### Changed

- Updated to CalculiX 2.20
- Reworked the support of Quasi 2D-3D simulations: second order elements are now allowed and all combination of write/read data are allowed ([#92](https://github.com/precice/calculix-adapter/pull/92)).

### Fixed

- Fixed several bugs and improved a few error messages. ([#88](https://github.com/precice/calculix-adapter/pull/88), [#94](https://github.com/precice/calculix-adapter/pull/94), [#95](https://github.com/precice/calculix-adapter/pull/95), [#102](https://github.com/precice/calculix-adapter/pull/102), [#103](https://github.com/precice/calculix-adapter/pull/103), [#104](https://github.com/precice/calculix-adapter/pull/104)).

## [v2.19.0] 2022-02-08

### Added

- Added Debian packages for Ubuntu 18.04 and 20.04 ([#72](https://github.com/precice/calculix-adapter/pull/72)).
- Added the flag `-USE_MT` to the Makefile, enabling multi-threading in SPOOLES ([#34](https://github.com/precice/calculix-adapter/pull/34), implemented in [#78](https://github.com/precice/calculix-adapter/pull/78)).
- Added experimental support for PaStiX ([#71](https://github.com/precice/calculix-adapter/pull/71) [#80](https://github.com/precice/calculix-adapter/pull/80) [#81](https://github.com/precice/calculix-adapter/pull/81)).
- Added clang-format and formatted the code ([#83](https://github.com/precice/calculix-adapter/pull/83)).

### Changed

- Updates for CalculiX v2.19 ([#78](https://github.com/precice/calculix-adapter/pull/)).

### Fixed

- Fixed some memory leaks ([#76](https://github.com/precice/calculix-adapter/pull/76) [#77](https://github.com/precice/calculix-adapter/pull/77)).
