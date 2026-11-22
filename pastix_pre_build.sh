#!/bin/bash
set -eu

ccx_src=${1:-"$HOME/CalculiX/ccx_2.24/src"}
mv "$ccx_src/pastix.h" "$ccx_src/pastix_ccx.h"
sed -i 's/#include "pastix.h"/#include "pastix_ccx.h"/g' "$ccx_src"/*.c
