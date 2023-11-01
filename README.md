# CalculiX-preCICE adapter

The adapter was initially developed for conjugate heat transfer (CHT) simulations via preCICE by Lucia Cheung in the scope of her master’s thesis [[1]](https://www5.in.tum.de/pub/Cheung2016_Thesis.pdf) in cooperation with [SimScale](https://www.simscale.com/). For running the adapter for CHT simulations refer to this thesis. The adapter was extended to fluid-structure interaction by Alexander Rusch [[2]](https://www.gacm2017.uni-stuttgart.de/registration/Upload/ExtendedAbstracts/ExtendedAbstract_0138.pdf).

The latest version of the adapter is based on **CalculiX version 2.21.**
The CalculiX adapter has been tested with the **preCICE version 2.5.0.**
Legacy versions of the adapter for older versions of CalculiX are supported on various branches. Branches for CalculiX version older than 2.15 require **preCICE v1.x**, whereas newer versions rely on **preCICE v2.x**.
 
## Start here

Go to the [adapter documentation](https://precice.org/adapter-calculix-overview.html) 

## References

[1] Lucia Cheung Yau. Conjugate heat transfer with the multiphysics coupling library precice. Master’s thesis, Department of Informatics, Technical University of Munich, 2016.

[2] Benjamin Uekermann, Hans-Joachim Bungartz, Lucia Cheung Yau, Gerasimos Chourdakis and Alexander Rusch. Official preCICE Adapters for Standard Open-Source Solvers. In Proceedings of the _7th GACM Colloquium on Computational Mechanics for Young Scientists from Academia_, 2017. 

## License

This project contains modified CalculiX files (licensed under [GPLv2 or later](http://www.dhondt.de/gpl-2.0.txt)) and additional adapter-specific files ([GPLv3](./LICENSE)). As a whole, this project is licensed under GPLv3.
