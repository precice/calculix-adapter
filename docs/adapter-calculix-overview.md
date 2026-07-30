---
title: The CalculiX adapter
permalink: adapter-calculix-overview.html
aliases:
  - /adapter-calculix-overview.html
  - /adapter-calculix.html
redirect_from: adapter-calculix.html
keywords: adapter, calculix
summary: "The CalculiX adapter can be used to couple CalculiX to CFD solvers for FSI or CHT application or even to couple CalculiX to itself."
---

## Start here

1. [Get CalculiX and the dependencies](adapter-calculix-get-calculix.html)
2. [Build the Adapter](adapter-calculix-get-adapter.html)
3. [Configure and run simulations](adapter-calculix-config.html)
4. Follow a [tutorial](https://precice.org/tutorials.html).

Are you encountering an unexpected error? Have a look at our [Troubleshooting](adapter-calculix-troubleshooting.html) page.

## Versions

The latest supported CalculiX version is {{site.calculix_version}}. If you already have a copy of the adapter, check the [adapter README](https://github.com/precice/calculix-adapter/blob/master/README.md) for the CalculiX version it was made for.

The adapter has a versioning scheme inherited from CalculiX: It is of the form `CCX_MAJOR.CCX_MINOR.ADAPTER_PATCH`. For instance, the release 2.20.0 modifies the source code of CalculiX 2.20. Further adapter releases for the same CalculiX version increase the `ADAPTER_PATCH` (e.g., 2.20.1), independent of whether it includes bug fixes, new features, or compatibility with a different preCICE version.

Compatibility with preCICE:

- preCICE v3 is supported since the adapter release `v2.20.1`
- preCICE v2 was supported in the branches `v2.16`, `v2.17` and the releases `v2.19.0` and `v2.20.0`
- preCICE v1 was supported in the branches `v2.10`, `v2.12`, `v2.13`, `v2.15`

## History

The adapter was initially developed for conjugate heat transfer (CHT) simulations via preCICE by Lucia Cheung in the scope of her master’s thesis[^1], in cooperation with [SimScale](https://www.simscale.com/). For running the adapter for CHT simulations refer to this thesis. The adapter was extended to fluid-structure interaction by Alexander Rusch[^2].

## References

[^1]: Lucia Cheung Yau. [Conjugate heat transfer with the multiphysics coupling library preCICE](https://mediatum.ub.tum.de/1461907). Master’s thesis, Department of Informatics, Technical University of Munich, 2016.

[^2]: Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. [Official preCICE Adapters for Standard Open-Source Solvers](https://doi.org/10.18419/opus-9334). In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017.
