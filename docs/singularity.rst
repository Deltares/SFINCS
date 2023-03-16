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

Using HPC Tooling to set up a SFINCS modelqueue
=========================================================

Premises
--------

The last part of the last seection was fine-and-dandy, but it does not use one thing of the HPC. In this section we set up a numnber of models for running
in the HPC's queuing mechanism.

Setup a model directory
-----------------------

Inside our environment, we set up a directory with models::

  $ tree models
  models
  ├── Blizard-Hazzard01
  ├── Blizard-Hazzard02
  ├── Charleston-subgrid01
  └── Charleston-subgrid02

Of course, **models** coould be any name here. Inside the subdirectory, you should have the usual sfincs configuration files, for example::

  $ tree Charleston-subgrid/
  Charleston-subgrid01
  ├── sfincs.bnd
  ├── sfincs.bzs
  ├── sfincs.dep
  ├── sfincs.dis
  ├── sfincs.ind
  ├── sfincs.inp
  ├── sfincs.msk
  ├── sfincs.obs
  ├── sfincs.sbg
  ├── sfincs.src
  └── sfincs.weir

Setting up the HPC Control Scripts
----------------------------------

The HPC Control Scripts really are two simple bash scripts which make it easier to run all the models in our **model** directory, submitting a job to a
HPC queueu for each model. 

The first shell script is a script to run a job (called test-job.sh, but you may call it whatever you want)::

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
  
This script isn't meant to run on its own and we get many errors if you attempt so (it is not user friendly either, since if you forget one parameter, it 
still tries something!). It is really to use for formatting an HPC job correctly. The second job is meant to use this script and submit the jobs to the queue (called sge-loop.sh)::

  #!/usr/bin/bash

  DIR=$1

  if [ -z ${DIR} ]
  then
    echo "Usage: $0 <model directory>"
    exit 0
  fi

  # Cleaning

  find . -name "*.nc" -delete
  find . -name "*.log" -delete
  find . -name "SFINCStest-*" -delete

  # Running qsub tasks

  for MODEL in `\ls $DIR`
  do
    echo "Running model in $MODEL"
    qsub -N SFINCStest-${MODEL} test-job.sh ${DIR} ${MODEL}
  done

  exit

This is also not very hard to understand. After **# Cleaning** there are some Unix commands that deletes old NetCDF output and old logfiles. Then it simply submits jobs to the HPC using **test-job.sh** as a template.


Tutorial run
============

Setup
-----

We have set some models up as follows::

  noorduin@v-hydrax001 ~/development/model-test $ find sums/
  sums/
  sums/Charleston-subgrid01
  sums/Charleston-subgrid01/sfincs.dis
  sums/Charleston-subgrid01/sfincs.bzs
  sums/Charleston-subgrid01/sfincs.msk
  sums/Charleston-subgrid01/sfincs.dep
  sums/Charleston-subgrid01/sfincs.sbg
  sums/Charleston-subgrid01/sfincs.bnd
  sums/Charleston-subgrid01/sfincs.src
  sums/Charleston-subgrid01/sfincs.weir
  sums/Charleston-subgrid01/sfincs.inp
  sums/Charleston-subgrid01/sfincs.ind
  sums/Charleston-subgrid01/sfincs.obs
  sums/Charleston-subgrid02
  sums/Charleston-subgrid02/sfincs.bnd
  sums/Charleston-subgrid02/sfincs.inp
  sums/Charleston-subgrid02/sfincs.ind
  sums/Charleston-subgrid02/sfincs.src
  sums/Charleston-subgrid02/sfincs.obs
  sums/Charleston-subgrid02/sfincs.weir
  sums/Charleston-subgrid02/sfincs.dep
  sums/Charleston-subgrid02/sfincs.sbg
  sums/Charleston-subgrid02/sfincs.bzs
  sums/Charleston-subgrid02/sfincs.msk

In this case it is rather lame (two times the same model), but you get the drift. I have used different names so that it is clear that the name **models** can also be changed. After this copy the shell scripts to th parent directory, so that things look like this::

  noorduin@v-hydrax001 ~/development/model-test $ ls -l

  drwxrwxr-x 5 noorduin domain users        5 Mar 16 08:39 models
  -rw-rw-r-- 1 noorduin domain users 66707456 Mar 14 11:59 sfincs-cpu.sif
  -rwxrwxr-x 1 noorduin domain users      354 Mar 16 08:45 sge-loop.sh
  -rwxrwxr-x 1 noorduin domain users      350 Mar 15 10:29 test-job.sh
  
Of course the **sfincs-cpu.sif** singualarity image can be in another directory, but then you have to adjust the **test-job.sh** script, so that it can find the image.

Run
---

This is simple::

noorduin@v-hydrax001 ~/development/model-test $ ./sge-loop.sh sums
Running model in Charleston-subgrid01
Your job 847944 ("SFINCStest-Charleston-subgrid01") has been submitted
Running model in Charleston-subgrid02
Your job 847945 ("SFINCStest-Charleston-subgrid02") has been submitted

and then::

  noorduin@v-hydrax001 ~/development/model-test $ qstat
  job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID
  -----------------------------------------------------------------------------------------------------------------
  ...
   847944 0.05500 SFINCStest noorduin     r     03/16/2023 09:07:10 test-c7@v-mcs055.directory.int     1
   847945 0.05500 SFINCStest noorduin     r     03/16/2023 09:07:10 test-c7@v-mcs056.directory.int     1
  ...

After a while you should see the usual error and info logs::

  noorduin@v-hydrax001 ~/development/model-test $ ls -l SFINCStest-*
  -rw-r--r-- 1 noorduin domain users  0 Mar 16 09:07 SFINCStest-Charleston-subgrid01.e847944
  -rw-r--r-- 1 noorduin domain users 74 Mar 16 09:07 SFINCStest-Charleston-subgrid01.o847944
  -rw-r--r-- 1 noorduin domain users  0 Mar 16 09:07 SFINCStest-Charleston-subgrid02.e847945
  -rw-r--r-- 1 noorduin domain users 74 Mar 16 09:07 SFINCStest-Charleston-subgrid02.o847945

and (for example)::

  noorduin@v-hydrax001 ~/development/model-test $ cat SFINCStest-Charleston-subgrid01.o847944
  v-mcs055.directory.intra 
  Running model Charleston-subgrid01.
  finished 


Exploring logfiles and result datasets
--------------------------------------

The system is designed so that logfiles as well as the result (NetCDF) datasets are written to the directory of the configuration of a model. So, for example::

  noorduin@v-hydrax001 ~/development/model-test $ ls -l sums/Charleston-subgrid01
  total 23048
  -rw-rw-r-- 1 noorduin domain users        0 Mar 16 09:07 error.log
  -rw-rw-r-- 1 noorduin domain users     2985 Mar 16 09:07 info.log
  -rwxrwxr-x 1 noorduin domain users      207 Nov 25  2021 sfincs.bnd
  -rwxrwxr-x 1 noorduin domain users   152179 Oct 24 11:29 sfincs.bzs
  -rwxrwxr-x 1 noorduin domain users   479456 Oct 29  2021 sfincs.dep
  -rwxrwxr-x 1 noorduin domain users    47848 Oct 24 11:29 sfincs.dis
  -rw-rw-r-- 1 noorduin domain users    10424 Mar 16 09:07 sfincs_his.nc
  -rwxrwxr-x 1 noorduin domain users   479460 Oct 29  2021 sfincs.ind
  -rwxrwxr-x 1 noorduin domain users      775 Nov  7 16:17 sfincs.inp
  -rw-rw-r-- 1 noorduin domain users 41832724 Mar 16 09:07 sfincs_map.nc
  -rwxrwxr-x 1 noorduin domain users   119864 Oct 29  2021 sfincs.msk
  -rwxrwxr-x 1 noorduin domain users      131 Oct 24 11:25 sfincs.obs
  -rwxrwxr-x 1 noorduin domain users 16301516 Oct 26  2021 sfincs.sbg
  -rwxrwxr-x 1 noorduin domain users       20 Oct 24 11:29 sfincs.src
  -rwxrwxr-x 1 noorduin domain users    59531 Oct 28  2021 sfincs.weir

Of course **sfincs_his.nc** and **sfincs_map.nc** are the result datasets, and **info.log** should contain the usual information about this run::

  noorduin@v-hydrax001 ~/development/model-test $ cat sums/Charleston-subgrid01/info.log

   ----------- Welcome to SFINCS -----------

    @@@@@  @@@@@@@ @@ @@  @@   @@@@   @@@@@
   @@@ @@@ @@@@@@@ @@ @@@ @@ @@@@@@@ @@@ @@@  
   @@@     @@      @@ @@@ @@ @@   @@ @@@
    @@@@@  @@@@@@  @@ @@@@@@ @@       @@@@@
       @@@ @@      @@ @@ @@@ @@   @@     @@@
   @@@ @@@ @@      @@ @@  @@  @@@@@@ @@@ @@@
    @@@@@  @@      @@ @@   @   @@@@   @@@@@

  ...

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
     5% complete,    26.1 s remaining ...
    10% complete,    24.9 s remaining ...
    ... 
    95% complete,     1.5 s remaining ...
   100% complete,     0.0 s remaining ...

   ---------- Simulation finished ----------

   Total time             :     30.172
   Total simulation time  :     30.084
   Time in input          :      0.088
   Time in boundaries     :      1.640 (  5.5%)
   Time in momentum       :     20.108 ( 66.8%)
   Time in continuity     :      8.159 ( 27.1%)
   Time in output         :      0.153 (  0.5%)

   Average time step (s)  :      6.878

   ---------- Closing off SFINCS -----------












