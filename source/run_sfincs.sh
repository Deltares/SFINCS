#!/bin/sh

LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH

# nvaccelinfo -v

export NVCOMPILER_ACC_TIME=0
export NVCOMPILER_ACC_SYNCHRONOUS=0

/usr/local/bin/sfincs/nvfortran_gpu_ccall/bin/sfincs
