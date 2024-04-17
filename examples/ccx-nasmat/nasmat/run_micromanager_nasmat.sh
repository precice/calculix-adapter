MKL_FOLDER=/opt/intel/oneapi/mkl/2023.1.0/lib/intel64
PETSC_DIR=/usr/lib/petsc
export LD_PRELOAD=$MKL_FOLDER/libmkl_def.so.2:$MKL_FOLDER/libmkl_avx2.so.2:$MKL_FOLDER/libmkl_core.so:$MKL_FOLDER/libmkl_intel_lp64.so:$MKL_FOLDER/libmkl_intel_thread.so:/usr/lib/x86_64-linux-gnu/libomp5.so:$PETSC_DIR/lib/libpetsc_real.so:$LD_PRELOAD

python3 run_micro_manager.py --config micro-manager-config.json
