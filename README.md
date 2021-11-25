# CalculiX-preCICE adapter

The adapter was initially developed for conjugate heat transfer (CHT) simulations via preCICE by Lucia Cheung in the scope of her master’s thesis [[1]](https://www5.in.tum.de/pub/Cheung2016_Thesis.pdf) in cooperation with [SimScale](https://www.simscale.com/). For running the adapter for CHT simulations refer to this thesis. The adapter was extended to fluid-structure interaction by Alexander Rusch [[2]](https://www.gacm2017.uni-stuttgart.de/registration/Upload/ExtendedAbstracts/ExtendedAbstract_0138.pdf).

This adapter was developed for **CalculiX version 2.16.** Other versions may be compatible, yet they have not been tested. Please let us know if you want to use a different version.

Adapters for other versions of CalculiX and preCICE are available in various branches. Branches compatible with **preCICE v2.x:**
 - master
 - v2.15_preCICE2.x
 
 All other branches are compatible with **preCICE v1.x**. 
 
## Start here

Go to the [adapter documentation](https://precice.org/adapter-calculix-overview.html) 

## References

[1] Lucia Cheung Yau. Conjugate heat transfer with the multiphysics coupling library precice. Master’s thesis, Department of Informatics, Technical University of Munich, 2016.

[2] Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. Official preCICE Adapters for Standard Open-Source Solvers. In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017. 

## License compatibility

The calculix-adapter is licensed under the [GPLv3](./LICENSE).
This is compatible with the GPLv2 License Calculix due to following statement:

From the official [license published on dhondt.de](http://www.dhondt.de/gpl-2.0.txt)

>    This program is free software; you can redistribute it and/or modify
>    it under the terms of the GNU General Public License as published by
>    the Free Software Foundation; either version 2 of the License, or
>    (at your option) any later version.
