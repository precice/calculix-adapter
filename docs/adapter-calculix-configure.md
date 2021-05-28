---
title: Configure the CalculiX adapter
permalink: adapter-calculix-config.html
keywords: adapter, calculix, configuration, config.yml
summary: "Write a config.yml, write a CalculiX case input file, and run an adapted CalculiX executable."
---

## Layout of the YAML configuration file

The layout of the YAML configuration file, which should be named `config.yml` (default name), is explained by means of an example for an FSI simulation:

```yaml
participants:

    Calculix:
        interfaces:
        - nodes-mesh: Calculix_Mesh
          patch: interface
          read-data: [Forces]
          write-data: [DisplacementDeltas]

precice-config-file: ../precice-config.xml
```

The adapter allows to use several participants in one simulation (e.g. several instances of Calculix if several solid objects are taken into account). The name of the participant "Calculix" must match the specification of the participant on the command line when running the executable of "CCX" with the adapter being used (this is described later). Also, the name must be the same as the one used in the preCICE configuration file *precice-config.xml*.  
One participant may have several FSI interfaces. Note that each interface specification starts with a dash.  
For FSI simulations the mesh type of an interface is always "nodes-mesh", i.e. the mesh is defined node-wise, not element-wise. The name of this mesh, "Calculix_Mesh", must match the mesh name given in the preCICE configuration file.  
For defining which nodes of the CalculiX domain belong to the FSI interface, a node set needs to be defined in the CalculiX input files. The name of this node set must match the name of the patch (here: "interface").  
In the current FSI example, the adapter reads forces from preCICE and feeds displacement deltas (not absolute displacements, but the change of the displacements relative to the last time step) to preCICE. This is defined with the keywords "read-data" and "write-data", respectively. The names (here: "Forces" and "DisplacementDeltas") again need to match the specifications in the preCICE configuration file. In the current example, the coupled fluid solver expects displacement deltas instead of displacements. However, the adapter is capable of writing either type. Just use "write-data: [Displacements]" for absolute displacements rather than relative changes being transferred in each time step. Valid `readData` keywords in CalculiX are:

```text
* Forces
* Displacements
* Temperature
* Heat-Flux
* Sink-Temperature
* Heat-Transfer-Coefficient
```

 Valid `writeData` keywords are:

```text
* Forces
* Displacements
* DisplacementDeltas
* Temperature
* Heat-Flux
* Sink-Temperature
* Heat-Transfer-Coefficient
```

From CalculiX version 2.15, additional `writeData` keywords are available:

```text
* Positions
* Velocities
```

Note that the square brackets imply that several read- and write-data types can be used on a single interface. This is not needed for FSI simulations (but for CHT simulations). Lastly, the "precice-config-file" needs to be identified including its location. In this example, the file is called *precice-config.xml* and is located one directory above the folder, in which the YAML configuration file lies.

## CalculiX case input file

An exemplary CalculiX case input file may look like the following:

```text
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

More input files that you may find in the CalculiX tutorial cases:

* `<name>.inp`: The main case configuration file. Through this, several other files are included.
* `<name>.msh`: The mesh file.
* `<name>.flm`: Films
* `<name>.nam`: Names, e.g. indices of boundary nodes
* `<name>.sur`: Surfaces
* `<name>.dfl`: DFlux

## Running the adapted calculiX executable

Running the adapted executable is pretty similar to running the original CalculiX CCX solver. The syntax is as follows:

```bash
ccx_preCICE -i [CalculiX input file] -precice-participant [participant name]
```

For example:

```bash
ccx_preCICE -i flap -precice-participant Calculix
```

The input file for this example would be *flap.inp*. Note that the suffix ".inp" needs to be omitted on the command line. The flag "-precice-participant" triggers the usage of the preCICE adapter. If the flag is not used, the original unmodified solver of CCX is executed. Therefore, the new executable "ccx_preCICE" can be used both for coupled preCICE simulations and CalculiX-only runs. Note that as mentioned above, the participant name used on the command line must match the name given in the YAML configuration file and the preCICE configuration file.

### Supported elements

The preCICE CalculiX adapter supports solid and shell elements. It can been used with both linear and quadratic tetrahedral (C3D4 and C3D10) and hexahedral (C3D8 and [C3D20](http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node29.html)) elements. For shell elements, currently S3 and S6 tetrahedral elements are supported. There is a restriction when using nearest-projection mapping that you have to use tetrahedral elements. If a quasi 2D-3D case is set up (single element in out-of-place direction) then only linear elements are supported.

### Nearest-projection mapping

In order to use nearest-projection mapping, a few additional changes are required. The first is that the interface surface file (.sur) must be added to the Calculix input file. An example of the addition to the input file is given below

```text
*INCLUDE, INPUT=all.msh
*INCLUDE, INPUT=fix1.nam
*INCLUDE, INPUT=fix2.nam
*INCLUDE, INPUT=fix3.nam
*INCLUDE, INPUT=interface.nam
*INCLUDE, INPUT=interface.sur
*MATERIAL, Name=EL
```

This surface file is generated during the mesh generation process. The second addition is to the config.yml. In order for the adapter to know that the surface mesh must be read, the line

```yaml
- nodes-mesh
```

must be changed to

```yaml
- nodes-mesh-with-connectivity
```

Note that an error will only occur if nodes-mesh-with-connectivity is specified without a .sur file. The calculix-adapter with nearest-projection mapping only supports tetrahedral elements (C3D4 and C3D10) as preCICE only works with surface triangles for nearest-projection mapping.
