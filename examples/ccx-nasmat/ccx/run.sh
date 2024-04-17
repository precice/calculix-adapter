#!/bin/sh
set -e -u

export OMP_NUM_THREADS=1
export CCX_NPROC_EQUATION_SOLVER=1
../../../bin/ccx_preCICE -i one_element -precice-participant macro_ccx
