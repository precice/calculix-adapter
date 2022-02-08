# NAME
ccx_preCICE - Adapter for using CalculiX with the preCICE library.

# SYNOPSIS
**ccx_preCICE** [*OPTION*] *files*

# DESCRIPTION

This page briefly summarizes how to use CalculiX with or without preCICE. For more detailed instructions, read the CalculiX documention on http://www.dhondt.de/ccx_2.17.pdf and preCICE documentation on https://precice.org/docs.html.
Omitting the -precice-participant option starts a stand-alone CalculiX simulation.


# OPTIONS

**-v, --version, -h, --help** 

: Prints current version of CalculiX. 

**-i < jobname >**

: Starts the job defined in <jobname>.inp. Be sure to omit the extension .inp in the job name.

**-precice-participant < participant >**

: Runs a preCICE simulation with name < participant >. This name must be the same as the one used in the adapter configuration file as well as in the preCICE configuration file.

# EXAMPLES

**ccx_preCICE -i flap -precice-participant CalculiX** : Runs the coupled simulation flap.inp as participant "CalculiX".

**ccx_preCICE -i flap** : Runs a CalculiX-only simulation with flap.inp as input file.

# REPORTING BUGS

: Report issues on the Github repository: https://github.com/precice/calculix-adapter.
