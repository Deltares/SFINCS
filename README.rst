=========================================
SFINCS: Super-Fast INundation of CoastS 
=========================================

|website| |docs_latest| |docs_stable| |license| |doi| |journal|

What is SFINCS?
-------------------------------------------------------
SFINCS is Deltares' new open-source, open-access reduced-complexity model designed for super-fast modelling of compound flooding events in a dynamic way!

Deltares is a not-for-profit independent institute for applied research in the field of water, subsurface and infrastructure. 

Why SFINCS?
-------------------------------------------------------
Compound flooding during extreme events can result in tremendous amounts of property damage and loss of life. Early warning systems and multi-hazard risk analysis can reduce these impacts. However, traditional approaches either do not involve relevant physics or are too computationally expensive to do so for large stretches of coastline. The SFINCS model is a new reduced-complexity engine recently developed by Deltares, that is capable of simulating compound flooding including a high computational efficiency balanced with good accuracy.

Where do I find more information about SFINCS?
-------------------------------------------------------
For general information see: https://www.deltares.nl/en/software/sfincs/

Find the user manual and more information on: https://sfincs.readthedocs.io/en/latest/

How do I get SFINCS?
-------------------------------------------------------
Download the latest windows executable here: https://download.deltares.nl/sfincs

Get the Docker of version of SFINCS to run on Mac, Linux or HPC here: https://hub.docker.com/r/deltares/sfincs-cpu

How to cite?
-------------------------------------------------------
To reference the software please use the the DOI provided in the SFINCS badge that points to the latest release: |doi|

The following paper presents the introduction of SFINCS:

   Leijnse, T., van Ormondt, M., Nederhoff, K., & van Dongeren, A. (2021). Modeling compound flooding in coastal systems using a computationally efficient reduced-physics solver: Including fluvial, pluvial, tidal, wind-      and wave-driven processes. Coastal Engineering, 163, 103796. https://doi.org/10.1016/j.coastaleng.2020.103796

How to contribute?
-------------------------------------------------------
If you find any issues in the code or documentation feel free to leave an issue on the `github issue tracker. <https://github.com/Deltares/SFINCS/issues>`_
You can find information about how to contribute to the SFINCS model at our `contributing page. <https://sfincs.readthedocs.io/en/latest/example.html#contributing>`_

SFINCS seeks active contribution from the hydro modelling community, so feel free to add something to our docs or model code, or reach out to 'sfincs@deltares.nl'!

How to ask questions to other users?
-------------------------------------------------------
For the SFINCS model itself, questions to fellow users can be asked in the SFINCS `github discussion board section. <https://github.com/Deltares/SFINCS/discussions>`_ 

Questions to fellow users regarding model building with the Python toolbox HydroMT-SFINCS can be asked in the HydroMT-SFINCS `github discussion board section. <https://github.com/Deltares/hydromt_sfincs/discussions>`_ 

Licensing conditions?
-------------------------------------------------------
This source code is available open source under GNU GPL-3.0 license, see attached.
This means (summarised, for details always read license):
	- Permissions: Commercial use, Distribution, Modification, Patent use, Private use
	- Conditions: Disclose source, License and copyright notice, Same license (so no distributing closed source versions), State changes
	- Limitations: No liability, no warranty
	
The SFINCS pre-compiled executables on Windows and Docker(linux) are available under the conditions as indicated when downloading SFINCS from https://download.deltares.nl/sfincs.

-For SFINCS versions v2.1.1 or OLDER-
SFINCS is available under Deltares free-trial license without monetary fee, meaning (summarised, for details always read license):
	- Permissions: No commercial use, Research publishing allowed, no distribution, no modification
	- Conditions: Disclose source, License and copyright notice, Same license (so no distribution of closed source versions), no changes
	- Limitations: No liability, no warranty
	- Support: The Software is provided to Licensee as is without any support from Deltares

-For SFINCS versions of 2025 or later-
SFINCS is available under Deltares Freeware license without monetary fee, meaning (summarised, for details always read license):
	- Permissions: Commercial use allowed, Research publishing allowed, no distribution, no modification
    	- IP, ownership: Deltares owns the copyright to SFINCS Software and the Deltares Software Documentation
	- Conditions: Disclose source, License and copyright notice, Same license (so no distribution of closed source versions), no changes
	- Limitations: No liability, no warranty
	- Support: Deltares will only provide any maintenance and/or support for SFINCS Software as part of separate project and/or a maintenance Agreement as agreed with Deltares

.. figure:: https://user-images.githubusercontent.com/28528822/200898347-d4016571-f3c7-4257-b59c-86aa1e97a699.png
   :width: 800px
   :align: center   
   
.. |website| image:: https://github.com/Deltares/SFINCS/blob/main/docs/figures/Deltares_logo_D-blauw_RGB.svg
    :target: https://www.deltares.nl/en/software-and-data/products/sfincs
    :alt: Website
    :width: 80px

.. |docs_latest| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
    :target: https://sfincs.readthedocs.io/en/latest
    :alt: Latest developers docs

.. |docs_stable| image:: https://img.shields.io/badge/docs-stable-brightgreen.svg
    :target: https://sfincs.readthedocs.io/en/v2.0.3_cauberg_release/
    :alt: Stable docs last release

.. |license| image:: https://img.shields.io/github/license/Deltares/SFINCS
    :alt: License
    :target: https://github.com/Deltares/SFINCS/blob/main/LICENSE    

.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10118583.svg
    :alt: Zenodo
    :target: https://doi.org/10.5281/zenodo.10118583

.. |journal| image:: https://github.com/Deltares/SFINCS/blob/main/docs/figures/SFINCS_logo.svg
    :alt: Elsevier
    :target: https://doi.org/10.1016/j.coastaleng.2020.103796    
    :width: 30px
