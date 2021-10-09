#! /bin/bash

echo "Renaming includes of \"pastix.h\" to \"pastix_ccx.h\" to avoid conflict."
sed -i 's/#include "pastix.h"/ #include "pastix_ccx.h"/g' /usr/local/CalculiX/ccx_2.17/src/*.c
echo "Renaming CalculiX file pastix.h to pastix_ccx.h"
mv /usr/local/CalculiX/ccx_2.17/src/pastix.h /usr/local/CalculiX/ccx_2.17/src/pastix_ccx.h
