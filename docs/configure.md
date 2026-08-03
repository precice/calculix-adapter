---
title: Configure the CalculiX adapter
permalink: adapter-calculix-config.html
aliases:
  - /adapter-calculix-config.html
keywords: adapter, calculix, configuration, config.yml
summary: "Write a config.yml, write a CalculiX case input file, and run an adapted CalculiX executable."
---

## Running the adapted CalculiX executable

Running the adapted executable is pretty similar to running the original CalculiX. The syntax is as follows:

```bash
ccx_preCICE -i [CalculiX input file] -precice-participant [participant name]
```

For example:

```bash
ccx_preCICE -i flap -precice-participant Calculix
```

The input file for this example would be `flap.inp`. Note that the suffix `.inp` needs to be omitted on the command line. The flag `-precice-participant` triggers the usage of the preCICE adapter. If the flag is not used, the original unmodified solver of CalculiX is executed, allowing CalculiX-only runs. Note that as mentioned above, the participant name used on the command line must match the name given in the YAML configuration file and in the preCICE configuration file.

## Adapter configuration file

The adapter looks for a YAML-based configuration file named `config.yml`, which, for historical reasons, starts by defining a list of participants. For example, for an FSI simulation:

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

The name of the participant `Calculix` must match the command-line argument `-precice-participant` for `CCX_preCICE` (read below) and the one used in the preCICE configuration file `precice-config.xml`. One participant may have several coupling interfaces. Note that each interface specification starts with a dash (new YAML list entry).

Depending on the data you need to read and write, the interface should define a mesh of one of the following types:

- a `faces-mesh` (or `mesh` as a synonym), where the data points are centers of faces (computed by the adapter). An interface made of faces should be defined in the CalculiX case using the `*SURFACE` command.
- a `nodes-mesh`, where the data points are the nodal vertices. An interface made of nodes should define these nodes using `*NSET`.
- an `elements-mesh`, where the data points are the quadrature points of the elements of a mesh. The mesh should be defined by nodes using `*NEST`. **Note**: `elements-mesh` is still experimental.

Using the wrong family of mesh (e.g. reading forces on faces) throws an error. If you need both kinds of meshes, you should define more than one interface.

In FSI simulations, the mesh type for an interface is always `nodes-mesh`, as forces and displacement are defined on nodes. The name of this mesh, `Calculix_Mesh`, must match the mesh name given in the preCICE configuration file. In CHT simulations, `faces-mesh` is used.
For defining which nodes of the CalculiX domain belong to the FSI interface, a node set needs to be defined in the CalculiX input files. The name of this node set must match the name of the patch (in this example, `interface`).

For multiscale mechanics simulations, the mesh type is always `elements-mesh`. The stresses, strains, and the material stiffness are defined on the quadrature points.

In this FSI example, the adapter reads forces from preCICE and feeds displacement deltas (not absolute displacements, but the change of the displacements relative to the last time step) to preCICE. This is defined with the keywords `read-data` and `write-data`, respectively. The names (here: `Forces` and `DisplacementDeltas`) again need to match the specifications in the preCICE configuration file. Absolute displacements can be configured with `Displacements`. 

Valid `readData` keywords in CalculiX are (with the corresponding boundary types; see the CalculiX documentation):

- On `faces-mesh`:
  - `Pressure` (Use a `*DLOAD`)
  - `Heat-Flux` (Use a `*DFLUX`)
  - `Sink-Temperature` (Use `*FILM`)
  - `Heat-Transfer-Coefficient` (Use `*FILM`)
- On `nodes-mesh`:
  - `Forces` (Use a `*CLOAD`)
  - `Displacements` (Use `*BOUNDARY`)
  - `Temperature` (Use `*BOUNDARY`)

Valid `writeData` keywords are:

- On `faces-mesh`:
  - `Pressure`
  - `Heat-Flux`
  - `Sink-Temperature`
  - `Heat-Transfer-Coefficient`
- On `nodes-mesh`:
  - `Forces`
  - `Displacements`
  - `DisplacementDeltas`
  - `Temperature`
  - `Positions`
  - `Velocities`

Note that the square brackets imply that several read- and write-data types can be used on a single interface (YAML list). This is mainly useful in CHT Robin coupling simulations.

Lastly, `precice-config-file` is the path to the preCICE configuration file.

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

When using `faces-meshes`, instead of a node set (`*NSET`), a `*SURFACE` must be sent, defined by a list of elements and face numbers. Instead of starting with an `N`, the name must start with a `S`.

CalculiX CCX offers both a geometrically linear and a geometrically non-linear solver, and both are supported by the adapter. The keyword `NLGEOM` (as shown in the example) selects the geometrically non-linear solver. It is also automatically triggered if material non-linearities are included in the analysis. In case the keyword `NLGEOM` does not appear in the CalculiX case input file and the chosen materials are linear, the geometrically linear CalculiX solver is used. In any case, for FSI simulations, the keyword `DYNAMIC` (enabling a dynamic computation) must appear in the CalculiX input file.

More input files that you may find in the CalculiX tutorial cases:

- `<name>.inp`: The main case configuration file. Through this, several other files are included.
- `<name>.msh`: The mesh file.
- `<name>.flm`: Films
- `<name>.nam`: Names, e.g., indices of boundary nodes
- `<name>.sur`: Surfaces
- `<name>.dfl`: DFlux

## Notes on meshes

### Supported elements

The CalculiX adapter supports most elements when using `nodes-mesh`. It has been used with both linear and quadratic tetrahedral (`C3D4` and `C3D10`) and hexahedral (`C3D8`, `C3D8I`, and `C3D20`) elements. For nearest-projection mapping, mesh connectivity is only provided when using tetrahedral elements.

For `faces-mesh`, tetrahedral and hexahedral meshes are supported.

When using `elements-mesh`, linear tetrahedral (C3D4) and hexahedral (C3D8) elements are supported.

### Coupling to 2D simulations

The adapter supports quasi-2D simulations when the z-direction is ignored. If you set the preCICE interface dimension to 2, the adapter will map data from the CalculiX 3D simulation to 2D space and vice-versa. The 3D simulation should be made of solid elements (or shells) of unit thickness.

For `nodes-mesh`, when writing continuous fields (such as temperature and displacements), the adapter will send data that is averaged over thickness. For conservative data (such as forces), sums are computed. When reading forces, the load applied to a 2D point will be spread evenly between the 3D points sharing the same x and y coordinates.

For `faces-mesh`, the z-component is discarded.

### Nearest-projection mapping

In order to use nearest-projection mapping, a few additional changes are required. The first is that the interface surface file (`.sur`) must be added to the CalculiX input file. For example:

```text
*INCLUDE, INPUT=all.msh
*INCLUDE, INPUT=fix1.nam
*INCLUDE, INPUT=fix2.nam
*INCLUDE, INPUT=fix3.nam
*INCLUDE, INPUT=interface.nam
*INCLUDE, INPUT=interface.sur
*MATERIAL, Name=EL
```

This surface file is generated during the mesh generation process. The second addition is to the config.yml. In order for the adapter to know that the surface mesh must be read, the adapter configuration file also needs to be modified:

```diff
- - nodes-mesh:
+ - nodes-mesh-with-connectivity:
```

Note that an error will only occur if `nodes-mesh-with-connectivity` is specified without a `.sur` file. The CalculiX adapter with nearest-projection mapping only supports tetrahedral elements (`C3D4` and `C3D10`), as preCICE only works with surface triangles for nearest-projection mapping.

## Modal dynamic simulations

The adapter supports modal dynamic simulations. In this type of simulation, eigenmodes from a frequency analysis are used. Instead of solving the full system of equations, CalculiX solves the problem as a time-dependent linear combination of these eigenmodes; this reduces the number of degrees of freedom of the system. Therefore, the simulation is faster, although its accuracy is dependent on the nonlinearity of the response. This method is very attractive in cases in which the solid dynamics is linear since, if you need many runs of the simulation, the `*FREQUENCY*` extraction step needs to be run only once.

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

## Parallelization

CalculiX comes with OpenMP and the SPOOLES library, which itself can use OpenMP. The adapter also supports this, and parallel runs can be used in the same way as with the uncoupled version of CalculiX. You can specify the number of threads via the `OMP_NUM_THREADS` environment variable. For a finer configuration, look at the CalculiX documentation.
You can also try [GPU acceleration with PaStiX](adapter-calculix-pastix-build.html).

## Restarting

To restart a CalculiX simulation, we need to enable restart files (`<name>.rout`), run the first simulation, and then restart the second simulation from that using a modified input file and renaming the `<name>.rout` to `<name>.rin`. See an [example from the community](https://pawel-lojek.medium.com/resuming-fsi-simulations-with-openfoam-calculix-896088861ae).

{% note %}
This section might be incomplete or contain inaccuracies. Help improve this page: Click "Edit me" to draft your suggestions.
{% endnote %}

1. In the `<name>.inp`, modify the end time in the following section (if needed):

   ```text
   *DYNAMIC, ALPHA=0.0, DIRECT
   1.E-2, 0.1
   ```

   The first number specifies the time step size, while the second specifies the duration of the current `STEP`. When restarting with the same number of time steps per `STEP`, the second number should not be modified.
2. Under the section specifying the time step size and end time, enable writing restart files (in this case, for every step):

   ```text
   *RESTART,WRITE,FREQUENCY=1
   ```

   At the very end of the simulation, and after a normal exit, a file `<name>.rout` will be generated.
   Rename this file to `<name>.rin`.
3. To restart a simulation, remove the mesh, material, and `*INCLUDE` sections from the input (`.inp`) file, keep/adapt the `*STEP` section(s), and add the following line as the first line of the file:
   
   ```text
   *RESTART,READ
   ```
  
   or
  
   ```text
   *RESTART,READ,STEP=1
   ```
   
   For every new step of restarting the simulation, increase the respective number: when you restart again to go further beyond in time, set `STEP=2`.
4. Since all the rest of the configuration is included in the restart file, we need to remove the rest of the definitions. In the end, the input file should look like this:

   ```text
   *RESTART,READ,STEP=1
   *STEP, INC=1000000
   *DYNAMIC, ALPHA=0.0, DIRECT
   1.E-2, 0.1
   *RESTART,WRITE,FREQUENCY=1
   *END STEP
   ```

   Note that, in CalculiX configuration files, you can add comment lines with two asterisks: `** comment`.
5. Make the respective adjustments in the other participants as well. For example, see [notes on restarting OpenFOAM simulations](https://precice.org/adapter-openfoam-config.html#restarting-fsi-simulations). Search also the [preCICE forum](https://precice.discourse.group/search?q=restarting%20order%3Alatest) and help improve this documentation.
