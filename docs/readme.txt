To build the docs locally:

- Make a local checkout of https://github.com/Deltares/hydromt_sfincs

- Open anaconda prompt (or other python environment manager)

- cd to the right folder: 'cd d:\repos\SFINCS_docs\docs\'

- activate environment: 'conda activate hydromt-sfincs'

- add sphinxcontrib package: conda install sphinxcontrib

- do: 'make html'

> The HTML pages are in _build\html.