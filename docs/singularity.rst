.. _docker-desktop:

Introduction
============

This documentation descibes the following:

o howto build the sfincs-cpu container on (for example the v-hydrax001.directory.intra). 

o howto automate the build process in TeamCity.

o In the "Tutorial" chapter of this document, I show how to run a simple simulation on the cluster with singularity.


Prerequisites on cluster user account
=====================================

To build a singularity cluster as well as using singularity runtime, you need to add the following modules information
to your **/home/${USER}/.bashrc file::

  # .bashrc

  # User specific aliases and functions

  alias rm='rm -i'
  alias cp='cp -i'
  alias mv='mv -i'

  # Source global definitions
  if [ -f /etc/bashrc ]; then
        . /etc/bashrc
  fi

  # Modules
  module load singularity
  module load squashfs-tools

When you logout and login, this file should be activated, and you can do::

  willem@deltares-test ~ $ singularity --version
  singularity-ce version 3.8.0


Building the singularity SFINCS container image
===============================================

In the future, this subsection should be covered in a Teamcity buildstep, which hopefully will spit out the sincs-cpu.sif
Also, a normal user doesn't have permission to build a Singularity container on the Rekencluster. That's why we  use a local
installed Virtual Machine. In reality, you should not have to do this becuase it is in a buildstep.

Cloning the SFINCS git repository
---------------------------------

This can be done by the following::

  willem@deltares-test ~/development/github $ git clone https://github.com/Deltares/SFINCS.git
  Cloning into 'SFINCS'...
  remote: Enumerating objects: 1520, done.
  remote: Counting objects: 100% (213/213), done.
  remote: Compressing objects: 100% (153/153), done.
  remote: Total 1520 (delta 102), reused 162 (delta 60), pack-reused 1307
  Receiving objects: 100% (1520/1520), 61.85 MiB | 5.12 MiB/s, done.
  Resolving deltas: 100% (820/820), done.
  Updating files: 100% (462/462), done.


Building the sfincs-cpu.sif image
---------------------------------

First checkout the right branch. In the near future, this **sfincs-singularity** branch is merged with **main**, so this should be fine then. For now::

  willem@deltares-test ~/development/github $ cd SFINCS/
  willem@deltares-test ~/development/github/SFINCS (main)$ git checkout remotes/origin/sfincs-singularity
  willem@deltares-test ~/development/github/SFINCS ((HEAD detached at origin/sfincs-singularity))$ git switch -c sfincs-singularity
  Switched to a new branch 'sfincs-singularity'

For building the image::

  willem@deltares-test:~/development/github/SFINCS/source$ cd source
  willem@deltares-test:~/development/github/SFINCS/source$ singularity build --fakeroot "/home/willem/development/singularity/sfincs/sfincs-cpu.sif" Singularityfile-cpu.def

  ... 
  <lots of output>
  ...

  INFO:    Creating SIF file...
  INFO:    Build complete: /home/willem/development/singularity/sfincs/sfincs-cpu.sif

and in the end we see::

  willem@deltares-test:~/development/github/SFINCS/source$ ls -l /home/willem/development/singularity/sfincs/sfincs-cpu.sif
  -rwxr-xr-x 1 willem willem 66621440 mrt  3 12:09 /home/willem/development/singularity/sfincs/sfincs-cpu.sif


Introducing the run-sfincs.sh script
------------------------------------

To easy a sfincs run, a shell-script **run-sfincs.sh** is provided. The script works like this::

  $ ./run-sfincs.sh
  Usage: ./run-sfincs.sh <directory in container which contains input>
  For example: ./run-sfincs.sh /data if the Singularity container is started as: ./sfincs-cpu --bind <input>:/data


Running a tutorial model inside sfincs on the HPC-cluster
=========================================================

Setup a run directory
---------------------

I have created a tar-ball of the sfinc-cpu.sif and a sample-model and copied this to the v-hydrax001. After untarring::

  noorduin@v-hydrax001 ~ $ find sfincs/
  sfincs/
  sfincs/Charleston-subgrid
  sfincs/Charleston-subgrid/sfincs.bnd
  sfincs/Charleston-subgrid/sfincs.bzs
  sfincs/Charleston-subgrid/sfincs.dep
  sfincs/Charleston-subgrid/sfincs.dis
  sfincs/Charleston-subgrid/sfincs.ind
  sfincs/Charleston-subgrid/sfincs.inp
  sfincs/Charleston-subgrid/sfincs.msk
  sfincs/Charleston-subgrid/sfincs.obs
  sfincs/Charleston-subgrid/sfincs.sbg
  sfincs/Charleston-subgrid/sfincs.src
  sfincs/Charleston-subgrid/sfincs.weir
  sfincs/sfincs-cpu.sif

Now for a run of this model::

	noorduin@v-hydrax001 ~/sfincs $ singularity shell -B /home/noorduin/sfincs/Charleston-subgrid/:/data sfincs-cpu.sif
	Singularity> run-sfincs.sh /data/
	
	 ----------- Welcome to SFINCS -----------
	
	  @@@@@  @@@@@@@ @@ @@  @@   @@@@   @@@@@
	 @@@ @@@ @@@@@@@ @@ @@@ @@ @@@@@@@ @@@ @@@
	 @@@     @@      @@ @@@ @@ @@   @@ @@@
	  @@@@@  @@@@@@  @@ @@@@@@ @@       @@@@@
	     @@@ @@      @@ @@ @@@ @@   @@     @@@
	 @@@ @@@ @@      @@ @@  @@  @@@@@@ @@@ @@@
	  @@@@@  @@      @@ @@   @   @@@@   @@@@@
	
	              ..............
	          ......:@@@@@@@@:......
	       ..::::..@@........@@.:::::..
	     ..:::::..@@..::..::..@@.::::::..
	    .::::::..@@............@@.:::::::.
	   .::::::..@@..............@@.:::::::.
	  .::::::::..@@............@@..::::::::.
	 .:::::::::...@@.@..@@..@.@@..::::::::::.
	 .:::::::::...:@@@..@@..@@@:..:::::::::..
	 ............@@.@@..@@..@@.@@............
	 ^^^~~^^~~^^@@..............@@^^^~^^^~~^^
	 .::::::::::@@..............@@.:::::::::.
	  .......:.@@.....@.....@....@@.:.......
	   .::....@@......@.@@@.@....@@.....::.
	    .:::~@@.:...:.@@...@@.:.:.@@~::::.
	     .::~@@@@@@@@@@.....@@@@@@@@@~::.
	       ..:~~~~~~~:.......:~~~~~~~:..
	          ......................
	              ..............
	
	 -----------------------------------------
	
	 Build-Revision: $Rev: v2.0.3-beta$
	 Build-Date:     $Date: 2023-02-24$
	
	 Reading input file ...
	 Info : Running SFINCS in subgrid mode ...
	 Reading meteo data ...
	 Info : Preparing SFINCS grid on regular mesh ...
	 Reading sfincs.ind ...
	 Reading sfincs.msk ...
	 Number of active z points    :       119864
	 Number of active u/v points  :       238715
	 Reading sfincs.sbg ...
	 Reading water level boundaries ...
	 Reading observation points ...
	 Initializing output ...
	
	 ---------- Starting simulation ----------
	
	   0% complete,       - s remaining ...
	   5% complete,    32.4 s remaining ...
	  10% complete,    31.7 s remaining ...
	  15% complete,    30.7 s remaining ...
	  20% complete,    29.4 s remaining ...
	  25% complete,    27.6 s remaining ...
	  30% complete,    25.6 s remaining ...
	  35% complete,    23.5 s remaining ...
	  40% complete,    21.7 s remaining ...
	  45% complete,    20.0 s remaining ...
	  50% complete,    18.9 s remaining ...
	  55% complete,    17.8 s remaining ...
	  60% complete,    16.4 s remaining ...
	  65% complete,    14.4 s remaining ...
	  70% complete,    12.3 s remaining ...
	  75% complete,    10.2 s remaining ...
	  80% complete,     8.1 s remaining ...
	  85% complete,     6.1 s remaining ...
	  90% complete,     4.1 s remaining ...
	  95% complete,     2.0 s remaining ...
	 100% complete,     0.0 s remaining ...
	
	 ---------- Simulation finished ----------
	
	 Total time             :     40.212
	 Total simulation time  :     40.122
	 Time in input          :      0.091
	 Time in boundaries     :      1.944 (  4.8%)
	 Time in momentum       :     27.052 ( 67.4%)
	 Time in continuity     :     10.951 ( 27.3%)
	 Time in output         :      0.146 (  0.4%)
	
	 Average time step (s)  :      6.878
	
	 ---------- Closing off SFINCS -----------
	Singularity> exit
	exit
	noorduin@v-hydrax001 ~/sfincs $ ls -l Charleston-subgrid/
	total 58128
	-rwxrwxr-x 1 noorduin domain users      207 Mar  3 12:29 sfincs.bnd
	-rwxrwxr-x 1 noorduin domain users   152179 Mar  3 12:29 sfincs.bzs
	-rwxrwxr-x 1 noorduin domain users   479456 Mar  3 12:29 sfincs.dep
	-rwxrwxr-x 1 noorduin domain users    47848 Mar  3 12:29 sfincs.dis
	-rw-rw-r-- 1 noorduin domain users    10424 Mar  3 12:42 sfincs_his.nc
	-rwxrwxr-x 1 noorduin domain users   479460 Mar  3 12:29 sfincs.ind
	-rwxrwxr-x 1 noorduin domain users      775 Mar  3 12:29 sfincs.inp
	-rw-rw-r-- 1 noorduin domain users 41832724 Mar  3 12:42 sfincs_map.nc
	-rwxrwxr-x 1 noorduin domain users   119864 Mar  3 12:29 sfincs.msk
	-rwxrwxr-x 1 noorduin domain users      131 Mar  3 12:29 sfincs.obs
	-rwxrwxr-x 1 noorduin domain users 16301516 Mar  3 12:29 sfincs.sbg
	-rwxrwxr-x 1 noorduin domain users       20 Mar  3 12:29 sfincs.src
	-rwxrwxr-x 1 noorduin domain users    59531 Mar  3 12:29 sfincs.weir
	noorduin@v-hydrax001 ~/sfincs $
	
Of course, the output are the **sfincs_his.nc** and the **sfincs_map.nc** files. If you run it like this, this output is 
written to the model directory and obtainable outside the **sfincs-cpu.sif** container.
