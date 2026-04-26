#!/bin/bash

MODELDIR=$1
MODEL=$2

#
# -- SGE options :
#

#$ -S /bin/bash
#$ -cwd
#$ -q test-c7
#$ -V

cd $SGE_O_WORKDIR

#
# -- the commands to be executed (programs to be run) :
#

echo $HOSTNAME 

# Modules
module load singularity

echo "Running model $MODEL."
singularity run -B ${MODELDIR}/${MODEL}:/data sfincs-cpu.sif
echo finished 

exit

