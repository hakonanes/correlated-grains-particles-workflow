# Correlated grain/particle analysis from directly combined EBSD/BSE data

Zenodo DOI for this repository coming.

This repository contains Jupyter notebook files, MATLAB scripts, and other files, apart from the raw BSE images and EBSD datasets, which are necessary to reproduce the results and figures in the paper *"Correlated subgrain and particle analysis of a recovered Al-Mn alloy by directly combining EBSD and backscatter electron imaging"*, which is published in *Materials Characterization* ([doi](https://doi.org/10.1016/j.matchar.2022.112228)):

```bibtex
@article{aanes2022correlated,
  author   = {Håkon W. Ånes and Antonius T. J. {van Helvoort} and Knut Marthinsen},
  title    = {Correlated subgrain and particle analysis of a recovered AlMn alloy by directly combining EBSD and backscatter electron imaging},
  doi      = {10.1016/j.matchar.2022.112228},
  issn     = {1044-5803},
  pages    = {112228},
  journal  = {Materials Characterization},
  keywords = {Electron backscatter diffraction, Texture analysis, Particle analysis, Image registration, Data fusion},
  year     = {2022},
}

```

The paper preprint is available on *arXiv* ([doi](https://doi.org/10.48550/arXiv.2205.05514)):

```bibtex
@article{anes2022correlated_arxiv,
  author  = {{\AA}nes, H{\aa}kon Wiik and van Helvoort, Antonius TJ and Marthinsen, Knut},
  title   = {Correlated subgrain and particle analysis of a recovered Al-Mn alloy by directly combining EBSD and backscatter electron imaging},
  doi     = {10.48550/arXiv.2205.05514},
  journal = {arXiv preprint arXiv:2205.05514},
  year    = {2022},
}
```

The raw EBSD and BSE data used in the paper is available on *Zenodo* ([doi](https://doi.org/10.5281/zenodo.6470217)):

```bibtex
@dataset{aanes2022correlated_data,
  author    = {Håkon Wiik Ånes and Antonius T. J. van Helvoort and Knut Marthinsen},
  title     = {{Electron backscatter diffraction data and backscatter electron images from a cold-rolled and recovered Al-Mn alloy}},
  doi       = {10.5281/zenodo.6470217},
  note      = {{The data was acquired while Håkon Wiik Ånes received financial support from the Norwegian University of Science and Technology (NTNU) through the NTNU Aluminium Product Innovation Centre (NAPIC).}},
  publisher = {Zenodo},
  year      = {2022},
}
```

The contents in this repository is licensed under the GPLv3+, since many of the softwares used have the same license.

## Contents

### notebooks

Static views of the notebooks are available via [nbviewer](https://nbviewer.org/github/hakonanes/correlated-grains-particles-workflow/tree/main/notebooks/).

Python packages used in the notebooks are listed in `requirements.txt` and can be installed into a virtual or conda environment:

```bash
pip install -r requirements.txt
```

* `ebsd1_dewrap.ipynb`: EBSD datasets acquired with a NORDIF detector are sometimes written with the last column of patterns as the first column. This notebook checks if this is the case, places the patterns correctly in the map, and writes them to an HDF5 file in the kikuchipy h5ebsd format.
* `ebsd2_preprocess.ipynb`: Generate indexing-independent views of an EBSD dataset (mean intensity map, virtual backscatter electron images, image quality map, and average neighbour dot product map) and calibrate the detector-sample geometry via projection center (PC) optimization with the [PyEBSDIndex Python package](https://github.com/USNavalResearchLaboratory/PyEBSDIndex). An average PC is used in dictionary indexing.
* `ebsd3_dictionary_indexing.ipynb`: Obtain crystal orientations from the EBSD patterns via dictionary indexing (DI) as implemented in kikuchipy. Requires an Al master pattern generated with EMsoft.
* `ebsd4_refinement.ipynb`: Refine crystal orientations obtained from DI.
* `generate_figures_in_paper.ipynb`: Generate final figures used in the paper.
* `image_registration.ipynb`: Manually identify control points in BSE images and EBSD intensity maps, perform image registration of BSE and EBSD datasets, and insert particles detected in BSE images and their sizes into EBSD datasets, making up the multimodal datasets.
* `particle_detection.ipynb`: Detect particles in BSE images and EBSD intensity maps.

Installation of packages has been tested to work on Linux (Ubuntu 22.04) and Windows 10. All notebooks have been tested to work on Linux, while the two indexing and refinement notebooks have also been tested to work on Windows 10.

### matlab_scripts

MATLAB packages used in the scripts are [MTEX](https://mtex-toolbox.github.io/) and [export_fig](https://mathworks.com/matlabcentral/fileexchange/23629-export_fig).

* `estimate_gnd.m`: Estimate geometrically necessary disloations prior to grain and texture analysis.
* `orientation_analysis.m`: Correlated analysis of subgrains and particles based on the multimodal dataset, specifically dispersoids at subgrain boundaries and subgrains at constituent particles.

### data

The control points and processed BSE images (reference images) and EBSD intensity maps (sensed images) used in image registration are included in this directory.
