!/bin/sh

#MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/man;
#MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/REDIST/cuda/12.4/compat;
#export MANPATH

#PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin:$PATH;
#PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/REDIST/cuda/12.4/compat:$PATH;
#export PATH

ldconfig /opt/nvidia/hpc_sdk/Linux_x86_64/24.5/REDIST/cuda/12.4/compat

export NVCOMPILER_ACC_TIME=1

cd /mnt/c/work/projects/sfincs/gpu_tests/parallel_benchmarks/02048

/usr/local/bin/sfincs/nvfortran_gpu_parallel_128_sync/bin/sfincs
