---
title: The CalculiX adapter
permalink: adapter-calculix-overview.html
redirect_from: adapter-calculix.html
keywords: adapter, calculix
summary: "The CalculiX adapter can be used to couple CalculiX to CFD solvers for FSI or CHT application or even to couple CalculiX to itself."
---

## Start here

1. [Get CalculiX and the dependencies](adapter-calculix-get-calculix.html)
2. [Build the Adapter](adapter-calculix-get-adapter.html)
3. [Configure and run simulations](adapter-calculix-config.html)
4. Follow a [tutorial](tutorials.html).

Are you encountering an unexpected error? Have a look at our [Troubleshooting](adapter-calculix-troubleshooting.html) page.

## Versions

The latest supported CalculiX version is {{site.calculix_version}}. If you already have a copy of the adapter, check the [adapter README](https://github.com/precice/calculix-adapter/blob/master/README.md) for the CalculiX version it was made for.

## History

The adapter was initially developed for conjugate heat transfer (CHT) simulations via preCICE by Lucia Cheung in the scope of her master’s thesis [[1]](https://mediatum.ub.tum.de/1461907) in cooperation with [SimScale](https://www.simscale.com/). For running the adapter for CHT simulations refer to this thesis. The adapter was extended to fluid-structure interaction by Alexander Rusch [[2]](https://doi.org/10.18419/opus-9334).

## References

[1] Lucia Cheung Yau. Conjugate heat transfer with the multiphysics coupling library precice. Master’s thesis, Department of Informatics, Technical University of Munich, 2016.

[2] Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. Official preCICE Adapters for Standard Open-Source Solvers. In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017.
