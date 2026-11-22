#!/bin/bash
set -eu

ccx_src=${1:-"$HOME/CalculiX/ccx_2.24/src"}
test -f "$ccx_src/pastix.h" || test -f "$ccx_src/pastix_ccx.h"
for source in "$ccx_src"/*.c; do
    if grep -Fq '#include "pastix.h"' "$source"; then
        sed -i.bak 's/#include "pastix.h"/#include "pastix_ccx.h"/g' "$source"
        rm "$source.bak"
    fi
done
if [ -f "$ccx_src/pastix.h" ]; then
    mv "$ccx_src/pastix.h" "$ccx_src/pastix_ccx.h"
fi
