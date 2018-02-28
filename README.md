# Building and Usage of the preCICE Adapter for CalculiX

This document describes how to build the preCICE adapter for CalculiX and how to use it for fluid-structure interaction (FSI) simulations. The adapter was initially developed for conjugate heat transfer (CHT) simulations via preCICE by Lucia Cheung in the scope of her master’s thesis [[1]](https://www5.in.tum.de/pub/Cheung2016_Thesis.pdf) in cooperation with [SimScale](https://www.simscale.com/). For running the adapter for CHT simulations refer to this thesis. The adapter was extended to fluid-structure interaction by Alexander Rusch.
This adapter was developed for CalculiX version 2.13. Other versions may be compatible, yet they have not been tested. Please let us know if you want to use a different version.

## Contents
<!-- toc orderedList:0 -->

- [Building and Usage of the preCICE Adapter for CalculiX](#building-and-usage-of-the-precice-adapter-for-calculix)
	- [Contents](#contents)
	- [Building the Adapter](#building-the-adapter)
        - [CalculiX](#calculix)
        - [preCICE](#precice)
        - [Adapter](#adapter)
	- [Running Simulations](#running-simulations)
		- [Layout of the YAML Configuration File](#layout-of-the-yaml-configuration-file)
        - [CalculiX Case Input File](#calculix-case-input-file)
        - [Running the Adapted CalculiX Executable](#running-the-adapted-calculix-executable)
	- [References](#references)

<!-- tocstop -->


## Building the Adapter
### CalculiX
Before installing the adapter, CalculiX itself must be downloaded and installed from http://www.dhondt.de. CalculiX consists of the solver, called "CCX" and a pre- and postprocessing software with graphical user interface "CGX". The installation procedure of CalculiX is described in the *README* files of the respective packages. Moreover, SPOOLES and ARPACK need to be installed for CalculiX. The procedure is also explained in the *README* files of CalculiX.

### preCICE
It is assumed that preCICE has been installed successfully beforehand. Concerning installation instructions for preCICE, have a look at the preCICE-wiki pages on GitHub: https://github.com/precice/precice/wiki/Building.

### Adapter
The adapter makes use of an additional input configuration file in YAML format. Therefore, the YAML parser "yaml-cpp" needs to be downloaded and installed from GitHub: https://github.com/jbeder/yaml-cpp. The building procedure is described in a *README* file included in the package.

The adapter source code can be downloaded from https://github.com/precice/calculix-adapter. The adapter is installed via a *Makefile*. It contains five variables that need to be adapted to your system:

 1. `CCX`: Location of the original CalculiX solver (CCX) source code ("src" directory)
 2. `SPOOLES`: Location of SPOOLES (top-level directory)
 3. `ARPACK`: Location of ARPACK (top-level directory)
 4. `PRECICE_ROOT`: Location of preCICE (top-level directory)
 5. `YAML`: Location of the YAML parser (configuration file reader) (top-level directory)

Furthermore, you might need to adapt the name of the ARPACK library `libarpack_INTEL.a` to your setting. 
Type "make" for building and "make clean" before a rebuild from the top-level directory of the adapter code, in which the *Makefile* is located. After building successfully, the executable "ccx_preCICE" is located in the "bin" folder.

## Running Simulations
### Layout of the YAML Configuration File
The layout of the YAML configuration file, which should be named *config.yml* (default name), is explained by means of an example for an FSI simulation:

```
participants:

    Calculix:
        interfaces:
        - nodes-mesh: Calculix_Mesh
          patch: interface
          read-data: [Forces0]
          write-data: [DisplacementDeltas0]

precice-config-file: ../precice-config.xml
```

The adapter allows to use several participants in one simulation (e.g. several instances of Calculix if several solid objects are taken into account). The name of the participant "Calculix" must match the specification of the participant on the command line when running the executable of "CCX" with the adapter being used (this is described later). Also, the name must be the same as the one used in the preCICE configuration file *precice-config.xml*.  
One participant may have several FSI interfaces. Note that each interface specification starts with a dash.  
For FSI simulations the mesh type of an interface is always "nodes-mesh", i.e. the mesh is defined node-wise, not element-wise. The name of this mesh, "Calculix_Mesh", must match the mesh name given in the preCICE configuration file.  
For defining which nodes of the CalculiX domain belong to the FSI interface, a node set needs to be defined in the CalculiX input files. The name of this node set must match the name of the patch (here: "interface").  
In the current FSI example, the adapter reads forces from preCICE and feeds displacement deltas (not absolute displacements, but the change of the displacements relative to the last time step) to preCICE. This is defined with the keywords "read-data" and "write-data", respectively. The names (here: "Forces0" and "DisplacementDeltas0") again need to match the specifications in the preCICE configuration file. Note that the "0" appending the data names is not optional. Since it is possible to define several interfaces, the data specifications need to be numbered starting from zero (e.g. Forces0, Forces1, Forces2, ...). In the current example, the coupled fluid solver expects displacement deltas instead of displacements. However, the adapter is capable of writing either type. Just use "write-data: [Displacements0]" for absolute displacements rather than relative changes being transferred in each time step. Note that the square brackets imply that several read- and write-data types can be used on a single interface. This is not needed for FSI simulations (but for CHT simulations). Lastly, the "precice-config-file" needs to be identified including its location. In this example, the file is called *precice-config.xml* and is located one directory above the folder, in which the YAML configuration file lies.

### CalculiX Case Input File
An exemplary CalculiX case input file may look like the following:

```
*INCLUDE, INPUT=all.msh
*INCLUDE, INPUT=fix1.nam
*INCLUDE, INPUT=fix2.nam
*INCLUDE, INPUT=fix3.nam
*INCLUDE, INPUT=interface.nam
*MATERIAL, Name=EL
*ELASTIC
 100000000, 0.3
*DENSITY
 10000.0
*SOLID SECTION, Elset=Eall, Material=EL
*STEP, NLGEOM, INC=1000000
*DYNAMIC
 0.01, 5.0
*BOUNDARY
 Nfix1, 3, 3, 0
 Nfix2, 1, 1, 0
 Nfix2, 3, 3, 0
 Nfix3, 1, 3, 0
*CLOAD
 Ninterface, 1, 0.0
 Ninterface, 2, 0.0
 Ninterface, 3, 0.0
*NODE FILE
 U
*EL FILE
 S, E
*END STEP
```

The adapter internally uses the CalculiX data format for point forces to apply the FSI forces at the coupling interface. This data structure is only initialized for those nodes, which are loaded at the beginning of a CalculiX analysis step via the input file. Thus, it is necessary to load all nodes of the node set, which defines the FSI interface in CalculiX (referring to the above example, the nodes of set "interface" (Note that in CalculiX a node set always begins with an "N" followed by the actual name of the set, which is here "interface".) are loaded via the "CLOAD" keyword.), in each spatial direction. However, the values of these initial forces can (and should) be chosen to zero, such that the simulation result is not affected.

CalculiX CCX offers both a geometrically linear as well as a geometrically non-linear solver. Both are coupled via the adapter. The keyword "NLGEOM" (as shown in the example) needs to be included in the CalculiX case input file in order to select the geometrically non-linear solver. It is also automatically triggered if material non-linearities are included in the analysis. In case the keyword "NLGEOM" does not appear in the CalculiX case input file and the chosen materials are linear, the geometrically linear CalculiX solver is used. In any case, for FSI simulations via preCICE the keyword "DYNAMIC" (enabling a dynamic computation) must appear in the CalculiX input file.

### Running the Adapted CalculiX Executable
Running the adapted executable is pretty similar to running the original CalculiX CCX solver. The syntax is as follows:

    ccx_preCICE -i [CalculiX input file] -precice-participant [participant name]

For example:

    ccx_preCICE -i flap -precice-participant Calculix

The input file for this example would be *flap.inp*. Note that the suffix ".inp" needs to be omitted on the command line. The flag "-precice-participant" triggers the usage of the preCICE adapter. If the flag is not used, the original unmodified solver of CCX is executed. Therefore, the new executable "ccx_preCICE" can be used both for coupled preCICE simulations and CalculiX-only runs. Note that as mentioned above, the participant name used on the command line must match the name given in the YAML configuration file and the preCICE configuration file.

## References
[1] Lucia Cheung Yau. Conjugate heat transfer with the multiphysics coupling library precice. Master’s thesis, Department of Informatics, Technical University of Munich, 2016.

[2] Alexander Rusch. Extending SU2 to fluid-structure interaction via preCICE. Bachelor's thesis, Munich School of Engineering, Technical University of Munich, 2016.
