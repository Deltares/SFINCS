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




For more information and docker-desktop itself, see https://www.docker.com/products/docker-desktop/.

Switching on necessary Windows capabilities
-------------------------------------------

Docker for Desktop is a quite complicated object (although the end result seems very obvious). Under the hood, a lot of virtual machinery
is going on. As a first, your laptops BIOS has to be configured so that the laptop is allowed to spawn virtual machines. Luckily, Deltares
laptops come with this preset already installed in the BIOS.

As a second, we have to switch on some windows capabilities. Although we can do this in the gui, the simplest is to do this via an elevated
Powershell prompt (so, running as Administrator). There, issue the following commands::

  PS C:\> Enable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform

(Answer N on the question whether to reboot)::

  PS C:> Enable-WindowsOptionalFeature -Online -FeatureName $("Microsoft-Hyper-V", "Containers") -All

(Answer Y on the question whether to reboot. This will reboot and you see that a couple of things get installed).


wsl 1
-----

Although it is not really necessary for GEOLibs, you want to include a Linux here, because Docker-Desktop expects it. We install an Ubuntu
in "Windows Subsystem for Linux 1" Here. From within the GUI, open en elevated Powershell. Then::

  PS C:\> wsl --list --online
  The following is a list of valid distributions that can be installed.
  The default distribution is denoted by '*'.
  Install using 'wsl --install -d <Distro>'.

    NAME               FRIENDLY NAME
  * Ubuntu             Ubuntu
    Debian             Debian GNU/Linux
    kali-linux         Kali Linux Rolling
    SLES-12            SUSE Linux Enterprise Server v12
    SLES-15            SUSE Linux Enterprise Server v15
    Ubuntu-18.04       Ubuntu 18.04 LTS
    Ubuntu-20.04       Ubuntu 20.04 LTS
    OracleLinux_8_5    Oracle Linux 8.5
    OracleLinux_7_9    Oracle Linux 7.9

  PS C:\> wsl --list --online
  The following is a list of valid distributions that can be installed.
  The default distribution is denoted by '*'.
  Install using 'wsl --install -d <Distro>'.

    NAME                                   FRIENDLY NAME
  * Ubuntu                                 Ubuntu
    Debian                                 Debian GNU/Linux
    kali-linux                             Kali Linux Rolling
    Ubuntu-18.04                           Ubuntu 18.04 LTS
    Ubuntu-20.04                           Ubuntu 20.04 LTS
    Ubuntu-22.04                           Ubuntu 22.04 LTS
    OracleLinux_8_5                        Oracle Linux 8.5
    OracleLinux_7_9                        Oracle Linux 7.9
    SUSE-Linux-Enterprise-Server-15-SP4    SUSE Linux Enterprise Server 15 SP4
    openSUSE-Leap-15.4                     openSUSE Leap 15.4
    openSUSE-Tumbleweed                    openSUSE Tumbleweed

  PS C:\> wsl --install -d Ubuntu
  Installing: Windows Subsystem for Linux
  Windows Subsystem for Linux has been installed.
  Installing: Windows Subsystem for Linux
  [==========================52,0%

After rebooting, the installation of Ubuntu in the Linux Subsystem for Windows will commence, and you see::

  Installing, this may take a few minutes...
  Please create a default UNIX user account. The username does not need to match your Windows username.
  For more information visit: https://aka.ms/wslusers
  Enter new UNIX username: deltares
  New password:
  Retype new password:
  passwd: password updated successfully
  Installation successful!
  To run a command as administrator (user "root"), use "sudo <command>".
  See "man sudo_root" for details.

  Welcome to Ubuntu 22.04.1 LTS (GNU/Linux 5.15.79.1-microsoft-standard-WSL2 x86_64)

  * Documentation:  https://help.ubuntu.com
  * Management:     https://landscape.canonical.com
  * Support:        https://ubuntu.com/advantage

  This message is shown once a day. To disable it please create the
  /home/deltares/.hushlogin file.

This is a ubuntu running inside Windows 10. Just like a normal Ubuntu::

  deltares@DESKTOP-ECMHLMF:~$ sudo su
  [sudo] password for deltares: ***

  root@DESKTOP-ECMHLMF:/home/deltares#
  ...


wsl 2
-----

The Docker Desktop is only running under wsl 2, so::

  PS C:\> wsl --update
  Checking for updates.
  The most recent version of Windows Subsystem for Linux is already installed.

  PS C:\> wsl --set-default-version 2

  PS C:\> wsl -l -v
    NAME      STATE           VERSION
  * Ubuntu    Stopped         2

  PS C:\Users\Willem> wsl -l -v
    NAME      STATE           VERSION
  * Ubuntu    Stopped         2


Installing Docker Desktop on Windows 10
---------------------------------------

This is a GUI Application, so go to the Windows 10 GUI and download and install https://www.docker.com/products/docker-desktop/. During
the installation, choose the installation option "Use WSL 2 instead of Hyper-V". Then let the installer do its work.

After login out and login in and starting docker desktop. You can do things like::

  PS C:\> docker ps
  CONTAINER ID   IMAGE     COMMAND   CREATED   STATUS    PORTS     NAMES
  PS C:\Users\Willem> docker images
  REPOSITORY   TAG       IMAGE ID   CREATED   SIZE
  PS C:\Users\Willem> docker search nginx
  NAME                                              DESCRIPTION                                     STARS     OFFICIAL   AUTOMATED
  nginx                                             Official build of Nginx.                        17941     [OK]
  linuxserver/nginx                                 An Nginx container, brought to you by LinuxSâ€¦   182
  bitnami/nginx                                     Bitnami nginx Docker Image                      150
 [OK]

Windows containers
------------------

The GEOLibs and the GEOApps (which is essentially GeoLibs bundleled with DStabilityConsole and DGeoFlowConsole) needs to run in windows
containers and that is an entirely different beast.

To achieve this, right-click the Docker icon in the task bar and choose "Switch to Windows Containers". Now we can do the following::

  PS C:\> hostname
  L02712              <== Outside the Container

  PS C:\> docker pull mcr.microsoft.com/windows/nanoserver:ltsc2019
  ltsc2019: Pulling from windows/nanoserver
  af0153d864f1: Pull complete
  Digest: sha256:fc2d54de31f170c0bef160137b4dc0a80c2105a218b248dc71c754e1fcabd14f
  Status: Downloaded newer image for mcr.microsoft.com/windows/nanoserver:ltsc2019
  mcr.microsoft.com/windows/nanoserver:ltsc2019

  PS C:\> docker images
  REPOSITORY                             TAG        IMAGE ID       CREATED       SIZE
  mcr.microsoft.com/windows/nanoserver   ltsc2019   00a00b91628a   2 weeks ago   258MB

  PS C:\> docker run -it mcr.microsoft.com/windows/nanoserver:ltsc2019 cmd
  Microsoft Windows [Version 10.0.17763.4010]
  (c) 2018 Microsoft Corporation. All rights reserved.

  C:\>hostname          <== Inside the Container
  ac77ce348bcc

  C:\>exit

  PS C:\> hostname
  L02712              <== Outside the Container
