% ccx_preCICE(1) 2.17
% Boris Martin
% Oct 2021

# NAME
ccx_preCICE - Adapter for using CalculiX with the preCICE library.

# SYNOPSIS
**ccx_preCICE** [*OPTION*] *files*

# DESCRIPTION

This page briefly summarizes how to use CalculiX with or without preCICE. For more detailed instructions, read the CalculiX documention on http://www.dhondt.de/ccx_2.17.pdf and preCICE documentation on https://precice.org/docs.html.
CalculiX is used alone when the -precice-participant option is not used.

## Using with preCICE

## Regular usage


# OPTIONS

**-v, --version, -h, --help** 

: Prints current version of CalculiX. 

**-i < jobname >**

: Starts the job defined in <jobname>.inp. Be sure to omit the extension .inp in the job name.

**-precice-participant < participant >**

: Runs a preCICE simulation with name < participant >. This name must be the one used in the YAML config file as well as in the preCICE configuration file.

# EXAMPLES

**ccx_preCICE -i flap -precice-participant Calculix** : Runs the coupled simulation flap.inp as participant "CalculiX".

**ccx_preCICE -i flap** : Runs a CalculiX-only simulation with flap.inp as input file.
