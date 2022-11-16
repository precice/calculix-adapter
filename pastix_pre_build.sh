#! /bin/bash

echo "Renaming includes of \"pastix.h\" to \"pastix_ccx.h\" to avoid conflict."
sed -i 's/#include "pastix.h"/ #include "pastix_ccx.h"/g' ~/CalculiX/ccx_2.20/src/*.c
echo "Renaming CalculiX file pastix.h to pastix_ccx.h"
mv ~/CalculiX/ccx_2.20/src/pastix.h ~/CalculiX/ccx_2.20/src/pastix_ccx.h
