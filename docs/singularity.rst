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

  noorduin@v-hydrax001 ~ $ singularity --version
  singularity-ce version 3.8.0


Building the singularity SFINCS container image
===============================================

In the future, this subsection should be covered in a Teamcity buildstep, which hopefully will spit out the sincs-cpu.sif

Cloning the SFINCS git repository
---------------------------------

This can be done by the following::

  noorduin@v-hydrax001 ~/development/github $ git clone https://github.com/Deltares/SFINCS.git
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

First checkout the right branch. In the near future, this **sfincs-singularity** branch is mergest with **main**, so this should be 
fine then. For now::

  noorduin@v-hydrax001 ~/development/github $ cd SFINCS/
  noorduin@v-hydrax001 ~/development/github/SFINCS (main)$ git checkout remotes/origin/sfincs-singularity
  noorduin@v-hydrax001 ~/development/github/SFINCS ((HEAD detached at origin/sfincs-singularity))$ git switch -c sfincs-singularity
  Switched to a new branch 'sfincs-singularity'

For building the image:






...


Tutorial: Running a simple SFINCS model
=======================================

...


