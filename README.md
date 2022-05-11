# Correlated grain/particle analysis from directly combined EBSD/BSE data

Zenodo DOI coming after first `release' of repository.

This repository contains Jupyter notebook files, MATLAB scripts, and other files, apart from the raw BSE images and EBSD datasets, which are necessary to reproduce the results and figures in the paper *"Correlated subgrain and particle analysis of a recovered Al-Mn alloy by directly combining EBSD and backscatter electron imaging"*, which was recently submitted to *Materials Characterization*. The paper preprint is available on [arXiv](). The raw data used in the paper, comprising three EBSD datasets and three sets of four BSE images, is available from Zenodo at [https://doi.org/10.5281/zenodo.6470217](https://doi.org/10.5281/zenodo.6470217).

## Contents

### notebooks

Static views of the notebooks are available via [`nbviewer`](https://nbviewer.org/github/hakonanes/correlated-grains-particles-workflow/notebooks).

Python packages used in the notebooks are listed in `requirements.txt` and can be installed into for example a conda environment:

```bash
conda create -n corr-env python=3.9
conda activate corr-env
pip install -r requirements.txt
```

* `ebsd1_dewrap.ipynb`: EBSD datasets acquired with a NORDIF detector are sometimes written with the last column of patterns as the first column. This notebook checks if this is the case, places the patterns correctly in the map, and writes them to an HDF5 file in the kikuchipy h5ebsd format.
* `ebsd2_preprocess.ipynb`: Generate indexing-independent views of an EBSD dataset (mean intensity map, virtual backscatter electron images, image quality map, and average neighbour dot product map) and calibrate the detector-sample geometry via projection center (PC) optimization with the [PyEBSDIndex Python package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex) (cubic materials only!). An average PC is used in dictionary indexing.
* `ebsd3_dictionary_indexing.ipynb`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in kikuchipy. Requires an Al master pattern generated with EMsoft.
* `ebsd4_refinement.ipynb`: Refine crystal orientations obtained from DI.
* `generate_figures_in_paper.ipynb`: Generate final figures used in the paper.
* `image_registration.ipynb`: Manually identify control points in BSE images and EBSD intensity maps, perform image registration of BSE and EBSD datasets, and insert particles detected in BSE images and their sizes into EBSD datasets, making up the multimodal datasets.
* `particle_detection.ipynb`: Detect particles in BSE images and EBSD intensity maps.

Installation of packages has been tested to work on Linux (Ubuntu 22.04) and Windows 10. All notebooks have been tested to work on Linux, while the two indexing and refinement notebooks have also been tested to work on Windows 10.

### matlab_scripts

MATLAB packages used in the scripts are [MTEX](https://mtex-toolbox.github.io/) and [export_fig](https://se.mathworks.com/matlabcentral/fileexchange/23629-export_fig).

* `estimate_gnd.m`: Estimate geometrically necessary disloations prior to grain and texture analysis.
* `orientation_analysis.m`: Correlated analysis of subgrains and particles based on the multimodal dataset, specifically dispersoids at subgrain boundaries and subgrains at constituent particles.

### data

The control points and processed BSE images (reference images) and EBSD intensity maps (sensed images) used in image registration are included in this directory. The raw EBSD datasets and BSE images are available from Zenodo at [https://doi.org/10.5281/zenodo.6470217](https://doi.org/10.5281/zenodo.6470217).
