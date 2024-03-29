#!/bin/bash
. /usr/share/Modules/init/bash

module load squashfs-tools
module load singularity

# Execute command singularity -d build --fakeroot <image>.sif recipe.def in order to get debug information
singularity build --fakeroot "sfincs_cpu_%singularityVersionNumberSFINCS%_lnx64_sif%build.counter%.sif" "Singularityfile-cpu.def"

tar -czvf sfincs_cpu_%singularityVersionNumberSFINCS%_lnx64_sif%build.counter%.tar.gz sfincs_cpu_%singularityVersionNumberSFINCS%_lnx64_sif%build.counter%.sif execute_singularity.sh readme.txt run_singularity.sh submit_singularity.sh

# Copy the artifact to network (for now to my homedirectory)
cp sfincs_cpu_%singularityVersionNumberSFINCS%_lnx64_sif%build.counter%.tar.gz /u/noorduin/from_tc/SFINCS/SFINCSset_lnx64_Singularity

