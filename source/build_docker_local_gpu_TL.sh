#!/bin/sh

cp -r /mnt/c/users/leijnse/repos/SFINCS/source /mnt/c/users/leijnse/tmp_build_sfincs

cd /mnt/c/users/leijnse/tmp_build_sfincs/source

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix
find . -name \*.sh|xargs dos2unix

echo "STARTING COMPILING SFINCS GPU $(date)"

docker build -f Dockerfile.gpu.25.5.ccall . -t leynse/sfincs-gpu:PR267_18032026 > build.log 2>&1

echo "FINISHED COMPILING SFINCS GPU $(date)"
