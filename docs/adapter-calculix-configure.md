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

The adapter allows us to use several participants in one simulation (e.g., several instances of Calculix if several solid objects are taken into account). The name of the participant `Calculix` must match the specification of the participant on the command line when running the executable of `CCX` with the adapter being used (this is described later). Also, the name must be the same as the one used in the preCICE configuration file `precice-config.xml`.  

One participant may have several coupling interfaces. Note that each interface specification starts with a dash.
Depending on the data you need to read and write, the interface should define either a `faces-mesh` (or `mesh` as a synonym) where the data points are centers of faces (computed by the adapter) or a mesh made of CalculiX vertices, with the keyword `nodes-mesh`. An interface made of faces should be defined in the CalculiX case using the `*SURFACE` command, whereas meshes with nodes should define these nodes using `*NSET`. Using the wrong family of mesh (e.g. reading forces on faces) throws an error. If you need both kinds of meshes, you should define more than one interface.

For FSI simulations the mesh type of an interface is always `nodes-mesh`, as forces and displacement are defined on nodes. The name of this mesh, `Calculix_Mesh`, must match the mesh name given in the preCICE configuration file. In CHT simulations, `faces-meshes` are usually chosen, as they are needed to apply heat fluxes or convective heat transfer.
For defining which nodes of the CalculiX domain belong to the FSI interface, a node set needs to be defined in the CalculiX input files. The name of this node set must match the name of the patch (here: "interface").  

In the current FSI example, the adapter reads forces from preCICE and feeds displacement deltas (not absolute displacements, but the change of the displacements relative to the last time step) to preCICE. This is defined with the keywords `read-data` and `write-data`, respectively. The names (here: `Forces` and `DisplacementDeltas`) again need to match the specifications in the preCICE configuration file. In the current example, the coupled fluid solver expects displacement deltas instead of displacements. However, the adapter is capable of writing either type. Just use `write-data: [Displacements]` for absolute displacements rather than relative changes being transferred in each time step. Valid `readData` keywords in CalculiX are:

On faces-mesh:

* Pressure (Use a `*DLOAD`)
* Heat-Flux (Use a `*DFLUX`)
* Sink-Temperature (Use `*FILM`)
* Heat-Transfer-Coefficient (Use `*FILM`)

On nodes-mesh:

* Forces (Use a `*CLOAD`)
* Displacements (Use `*BOUNDARY`)
* Temperature (Use `*BOUNDARY`)

Have a look at the CalculiX documentation for a detailed description of each of these commands. There is an [online (but outdated) version](https://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node1.html) and an [up-to-date PDF version](http://www.dhondt.de/ccx_2.19.pdf).

 Valid `writeData` keywords are:

On faces-mesh:

* Pressure
* Heat-Flux
* Sink-Temperature
* Heat-Transfer-Coefficient

On nodes-mesh:

* Forces
* Displacements
* DisplacementDeltas
* Temperature

From CalculiX version 2.15, additional `writeData` keywords are available:

```text
* Positions
* Velocities
```

Note that the square brackets imply that several read- and write-data types can be used on a single interface. This is not needed for FSI simulations (but for CHT simulations). Lastly, the `precice-config-file` needs to be identified including its location. In this example, the file is called `precice-config.xml` and is located one directory above the folder, in which the YAML configuration file lies.

## CalculiX case input file

CalculiX is designed to be compatible with the Abaqus file format. Here is an example of a CalculiX input file:

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

The adapter internally uses the CalculiX data format for point forces to apply the FSI forces at the coupling interface. This data structure is only initialized for those nodes, which are loaded at the beginning of a CalculiX analysis step via the input file. Thus, it is necessary to load all nodes of the node set, which defines the FSI interface in CalculiX (referring to the above example, the nodes of set `interface` (Note that in CalculiX a node set always begins with an `N` followed by the actual name of the set, which is here `interface`.) are loaded via the `CLOAD` keyword.), in each spatial direction. However, the values of these initial forces can (and should) be chosen to zero, such that the simulation result is not affected.

When using "faces-meshes", instead of a node set (`*NSET`), a `*SURFACE` must be sent, defined by a list of elements and face numbers. Instead of starting with a "N", the name must start with a "S".

CalculiX CCX offers both a geometrically linear as well as a geometrically non-linear solver. Both are coupled via the adapter. The keyword "NLGEOM" (as shown in the example) needs to be included in the CalculiX case input file in order to select the geometrically non-linear solver. It is also automatically triggered if material non-linearities are included in the analysis. In case the keyword "NLGEOM" does not appear in the CalculiX case input file and the chosen materials are linear, the geometrically linear CalculiX solver is used. In any case, for FSI simulations via preCICE the keyword "DYNAMIC" (enabling a dynamic computation) must appear in the CalculiX input file.

More input files that you may find in the CalculiX tutorial cases:

* `<name>.inp`: The main case configuration file. Through this, several other files are included.
* `<name>.msh`: The mesh file.
* `<name>.flm`: Films
* `<name>.nam`: Names, e.g. indices of boundary nodes
* `<name>.sur`: Surfaces
* `<name>.dfl`: DFlux

## Running the adapted CalculiX executable

Running the adapted executable is pretty similar to running the original CalculiX CCX solver. The syntax is as follows:

```bash
ccx_preCICE -i [CalculiX input file] -precice-participant [participant name]
```

For example:

```bash
ccx_preCICE -i flap -precice-participant Calculix
```

The input file for this example would be `flap.inp`. Note that the suffix `.inp` needs to be omitted on the command line. The flag `-precice-participant` triggers the usage of the preCICE adapter. If the flag is not used, the original unmodified solver of CalculiX is executed. Therefore, the new executable `ccx_preCICE` can be used both for coupled preCICE simulations and CalculiX-only runs. Note that as mentioned above, the participant name used on the command line must match the name given in the YAML configuration file and the preCICE configuration file.

### Supported elements

The preCICE CalculiX adapter should support most elements when using `nodes-meshes`. It has been used with both linear and quadratic tetrahedral (C3D4 and C3D10) and hexahedral (C3D8, C3D8I, and [C3D20](http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node29.html)) elements. There is however a restriction when using nearest-projection mapping: in that case, you have to use tetrahedral elements.

When using face meshes, only tetrahedra and hexaedra are supported.

### Coupling to 2D simulations

The adapter supports quasi 2D simulations when the z-direction is ignored. If you set the preCICE interface dimension to 2, the adapter will map data from the CalculiX 3D simulation to 2D space and vice-versa. The 3D simulation should be made of solid elements (or shells) of unit thickness.

#### Behavior with `nodes-mesh`

When writing continuous fields (such as temperature and displacements), the adapter will send data that is averaged over thickness. For conservative data (such as forces), sums are computed. When reading forces, the load applied to a 2D point will be spread evenly between the 3D points sharing the same x and y coordinates.

#### Behavior with `faces-mesh`

When using `faces-mesh`, the behavior in unchanged and the z-component is discarded.

### Nearest-projection mapping

In order to use nearest-projection mapping, a few additional changes are required. The first is that the interface surface file (`.sur`) must be added to the CalculiX input file. An example of the addition to the input file is given below

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

Note that an error will only occur if nodes-mesh-with-connectivity is specified without a `.sur` file. The CalculiX adapter with nearest-projection mapping only supports tetrahedral elements (C3D4 and C3D10) as preCICE only works with surface triangles for nearest-projection mapping.

### Modal dynamic simulations

The adapter supports modal dynamic simulations. In this type of simulation, eigenmodes from a frequency analysis are used. Instead of solving the full system of equations, CalculiX solves the problem as a time dependent linear combination of these eigenmodes; this reduces the number of degrees of freedom of the sytem. Therefore, the simulation is faster, although its accuracy is dependent on the nonlinearity of the response. This method is very attractive in cases in which the solid dynamics is linear since, if you need many runs of the simulation, the `*FREQUENCY*` extraction step needs to be run only once.

To run a case with a `*MODAL DYNAMIC` analysis, the adapter has a special requirement. The `*FREQUENCY` step and the `*MODAL DYNAMIC` step must be run in different input files, as it is a requirement of the adapter to extract the frequency and modal data with a dedicated `.inp` file and then run the modal dynamic analysis without adding the frequency extraction step. The first run creates a `.eig` file with the modal information needed for the second run.

To run the frequency analysis without preCICE, create an input file like the following example:

```text
*INCLUDE, INPUT=all.msh
*INCLUDE, INPUT=fix1.nam
*INCLUDE, INPUT=fix2.nam
*INCLUDE, INPUT=fix3.nam
*INCLUDE, INPUT=interface.nam
*INCLUDE, INPUT=interface.sur
**===============================================================
*MATERIAL, Name=EL
*ELASTIC
 100000000, 0.3
*DENSITY
 10000.0
*BOUNDARY
 Nfix1, 3, 3, 0
 Nfix2, 1, 1, 0
 Nfix2, 3, 3, 0
 Nfix3, 1, 3, 0
**===============================================================
*STEP
*FREQUENCY, STORAGE=YES
4
*CLOAD
 Ninterface, 1, 0.0
 Ninterface, 2, 0.0
 Ninterface, 3, 0.0
*NODE FILE
 U
*EL FILE
 S, E
*END STEP
**===============================================================
```

Run the frequency analysis with:

```bash
ccx_preCICE -i [CalculiX input file]
```

Now the modal information should be stored in the `.eig` file. For the regular coupled simulation, delete the `*FREQUENCY` step and add the `*MODAL DYNAMIC` step as shown in the following example:

```text
*INCLUDE, INPUT=all.msh
*INCLUDE, INPUT=fix1.nam
*INCLUDE, INPUT=fix2.nam
*INCLUDE, INPUT=fix3.nam
*INCLUDE, INPUT=interface.nam
*INCLUDE, INPUT=interface.sur
**===============================================================
*MATERIAL, Name=EL
*ELASTIC
 100000000, 0.3
*DENSITY
 10000.0
*BOUNDARY
 Nfix1, 3, 3, 0
 Nfix2, 1, 1, 0
 Nfix2, 3, 3, 0
 Nfix3, 1, 3, 0
**===============================================================
*STEP, INC=1000000000
*MODAL DYNAMIC
1.E-4, 1
*CLOAD
 Ninterface, 1, 0.0
 Ninterface, 2, 0.0
 Ninterface, 3, 0.0
*NODE FILE
 U
*EL FILE
 S, E
*END STEP
**===============================================================
```

Run the standard coupled simulation with:

```bash
ccx_preCICE -i [CalculiX input file] -precice-participant [CalculiX participant name]
```

Make sure to replace `[CalculiX input file]` and `[CalculiX participant name]` with the appropriate file and name for your case.

### Parallelization

CalculiX comes with OpenMP and the SPOOLES library which itself can use OpenMP. The adapter also supports this and parallel runs can be used in the same way as with the uncoupled version of CalculiX. You can specify the number of threads via the `OMP_NUM_THREADS` environment variable. For a finer configuration, look at the CalculiX documentation.
You can also try [GPU acceleration with PaStiX](adapter-calculix-pastix-build.html).
