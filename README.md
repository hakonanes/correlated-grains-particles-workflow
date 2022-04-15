# Correlated grain/particle analysis from combined EBSD/BSE data

Zenodo DOI coming after first `release' of repository.

This repository contains all the Jupyter notebook files, MATLAB scripts, and other files, apart from the raw BSE images and EBSD data sets, to reproduce the results and figures in the paper *"Correlated subgrain and particle analysis in a recovered Al-Mn alloy by directly combining EBSD and backscatter electron imaging"*, which was recently submitted to the journal ?. The paper preprint is available on [arXiv]().

## Contents

### notebooks

Static views of the notebooks are available via [`nbviewer`](https://nbviewer.org/github/hakonanes/correlated-grains-particles-workflow/notebooks).

Python packages used in the notebooks are listed in `requirements.txt` and can be installed into for example a conda environment:

```bash
conda create -n corr-env python=3.9
conda activate corr-env
pip install -r requirements.txt
```

* `ebsd1_dewrap.ipynb`: EBSD data sets acquired with a NORDIF detector are sometimes written with the last column of patterns as the first column. This notebook checks if this is the case, places the pattern correctly in the map, and writes them to an HDF5 file in the kikuchipy h5ebsd format.
* `ebsd2_preprocess.ipynb`: Generate indexing-independent views of an EBSD data set (mean intensity map, virtual backscatter electron images, image quality map, and average neighbour dot product map) and calibrate the detector-sample geometry via projection center (PC) optimization with the [PyEBSDIndex Python package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex) (cubic patterns only!). An average PC is used in dictionary indexing.
* `ebsd3_dictionary_indexing.ipynb`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in kikuchipy. Requires an Al master pattern generated with EMsoft.
* `ebsd4_refinement.ipynb`: Refine crystal orientations obtained from DI.
* `generate_figures_in_paper.ipynb`: Generate final figures used in the paper.
* `image_registration.ipynb`: Manually identify control points in BSE images and EBSD intensity maps, perform image registration of BSE and EBSD data sets, and insert particles detected in BSE images into EBSD data sets.
* `particle_detection.ipynb`: Detect particles in BSE images and EBSD intensity maps.

### matlab_scripts

MATLAB packages used in the scripts are [MTEX](https://mtex-toolbox.github.io/) and [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig).

* `estimate_gnd.m`: Estimate geometrically necessary disloations prior to grain and texture analysis.
* `orientation_analysis.m`: Correlated analysis of subgrains and particles in the (combined) multi-modal data set, specifically subgrains by constituent particles and dispersoids on subgrain boundaries.

### data

Data produced from the Jupyter notebooks and MATLAB scripts of the three data sets with EBSD data and a BSE image.
