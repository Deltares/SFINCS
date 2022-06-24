Setting up models
=======

Introduction 
----------------------

SFINCS models can be setup using simple ascii text and binary input files, which can be generated on whatever platform suiting you as a user best.
To make the setting up of basic models easier, making a SFINCS model is supported by the two easy setup tools called 'Delft Dashboard' and 'HydroMT'.
These options do not cover setting up subgrid feature type models yet, which is still an advanced user functionality.
Additionaly, a range of Matlab scripts exist for setting up SFINCS models using the Open Earth Tools for more user flexibility.
If you need a more tailor-made solution for setting up your SFINCS models get in touch with us!

Delft Dashboard 
----------------------

Delft dashboard (Van Ormondt et al. 2020 https://doi.org/10.2166/hydro.2020.092) is a quick set-up tool for hydrodynamic models that includes setting up SFINCS models.
The tool can be run on Matlab or as standalone executable and has all the basic functionality to setup your basic model using globally available DEMs.
For more information see: https://publicwiki.deltares.nl/display/DDB

HydroMT 
----------------------

HydroMT is a more recent addition to the tools available for setting up SFINCS models, and is a Python based alternative.
Besides globally available DEMs it can also retrieve spatially varying infiltration and manning roughness data based on landuse maps.
Also, it is possible to setup a offline coupled model together with the hydrological wflow model that will provide boundary conditions of river discharge.
For more information regarding the SFINCS plugin of HydroMT see: https://deltares.github.io/hydromt_sfincs/
For more information regarding HydroMT in general see: https://deltares.github.io/hydromt/
For more user flexibility, it is also possible now to access individual setup components to build your own model or forcing from scratch, see: https://deltares.github.io/hydromt_sfincs/latest/user_guide/sfincs.html
For an example of building a model from python see: https://deltares.github.io/hydromt_sfincs/latest/_examples/build_from_py.html

Open Earth Tools
----------------------

In the Open Earth Tools section of DDB, a number of Matlab scripts are included to have more flexibility in setting up your SFINCS models.
Throughout the User Manual examples of how to use these scripts are given in the code blocks '**Matlab example using OET**'.
For a general overview of possible input files and how to create them using Matlab scripts see: https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/DelftDashBoard/general/sfincs/fileio/Input_format_examples.m
For getting started with Open Earth Tools see: https://publicwiki.deltares.nl/display/OET/OpenEarth

Running SFINCS
=======

SFINCS can be run on multiple different platforms, both local and cloud based.
The most commonly used and most extensively tested method is running SFINCS on windows using a batch-file.

On windows (standard)
----------------------

The standard method to run SFINCS locally is on a windows machine using a batch-file.
This batch file you copy to the folder where your input files to be used by SFINCS are located.
The batch file simply calls the executable (add the right path to the folder where sfincs.exe is located) and the general output text file is written to a new text file called 'sfincs_log.txt', see below for an example.

**To get a copy of the SFINCS executable get in touch with us, this is currently only possible when signing a pre-release software agreement agreeing that the software is property of Deltares and cannot be used for commercial purposes, only for testing of the code untill SFINCS has been made completely open source!**

Using batch-file
^^^^^^^^^

**run.bat**

.. code-block:: text	
	
	make a text file called 'run.bat' and add here:
	
		call "c:\..\folder_where_exe_is_located\sfincs.exe">sfincs_log.txt	
	
On linux 
----------------------

With windows exe using Wine
^^^^^^^^^

Some tests have been performed succesfully to run the SFINCS windows executable on a desktop linux system using Wine (https://www.winehq.org/).
For more questions regarding this option get in touch.


Dedicated linux compiled version
^^^^^^^^^

We are still performing tests to support the option of a dedicated linux compiled version of SFINCS, for more questions regarding this option get in touch.


Using Docker
----------------------

For always using the last build version of SFINCS on Windows, Mac or a cloud based cluster a convenient solution is running a Docker container version of SFINCS.
This can be done on a local desktop or in a cloud based cluster supporting docker (or using singularity, see below).

**Note that one is only allowed to use the online Docker version of SFINCS ONLY after signing a pre-release software agreement agreeing that the software is property of Deltares and cannot be used for commercial purposes, only for testing of the code untill SFINCS has been made completely open source! Get in touch to arrange this license agreement.**

Local desktop version
^^^^^^^^^

After downloading Docker desktop for your operating system (https://www.docker.com/products/docker-desktop), you can run a model using:

**Example**

.. code-block:: text

	docker pull deltares/sfincs-cpu

	docker run -vC:/Users/../SFINCS:/data deltares/sfincs-cpu

	(here 'C:/Users/../SFINCS' is the folder where the SFINCS input files to be used are located)

Cloud based cluster
^^^^^^^^^

The same principle is also possible on a cloud based cluster that supports running docker containers

Using Singularity
----------------------

On cloud based clusters like Surfsara/Azure/Amazon that **supports singularity**, it is possible to run the Docker container version of SFINCS directly.
Depending on the application it could be wise to pull the docker container once and save as new image, after which this image can be run multiple times.
This prevents unnesissarily loading the Docker container every time a simulation is performed.

**Note that one is only allowed to use the online Docker version of SFINCS ONLY after signing a pre-release software agreement that the software is property of Deltares and cannot be used for commercial purposes, only for testing of the code untill SFINCS has been made completely open source! Get in touch to arrange this license agreement.**

**Example**

.. code-block:: text	
	
	Pulling and running the docker container immediately:
	
		singularity run -B$(pwd):/data --nv docker://deltares/sfincs-cpu

	
	First pulling the docker container and creating a singularity image, then running this image:
	
		singularity pull docker://deltares/sfincs-cpu sfincs-cpu.img

		singularity run -B$(pwd):/data sfincs-cpu.img
	
	
Contributing
=======

Documentation 
----------------------	

The docs code has been moved to https://github.com/Deltares/SFINCS as of 24-06-2022.
Get in touch if you have suggestions how to improve this manual!

Code 
----------------------	

The SFINCS code as of now is on a closed community bases, get in touch if you would like to join us in developing the SFINCS code!
