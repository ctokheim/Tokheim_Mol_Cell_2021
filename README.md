# Genomic Landscape of the Ubiquitin-Proteasome System in Human Cancers

The Ubiquitin-Proteasome System (UPS), with over 600 proteins in its pathway, is responsible for the degradation of over 80% of human proteins. While UPS-targeting cancer drugs are under rapid development, a clearer understanding of which UPS genes are implicated in human cancer development and progression with their downstream protein substrates is necessary. By leveraging multi-omics data across more than 9,000 human tumors and 33 cancer types, we found that approximately 19% of all cancer driver genes impact the function of the UPS and that transcription factors are important substrates. Moreover, we developed a deep learning model (deepDegron) to identify substrate-intrinsic mechanisms where degron loss leads to up-regulated protein abundance, such as GATA3 and PPM1D. The UPS driver genes are substantially correlated with immune-related biomarkers, including IFNG response and lymphocyte infiltration, and overall patient survival. This study underlies the important role of UPS dysregulation in human cancers and the potential therapeutic utilities of targeting the UPS.

This repository contains jupyter notebooks used for data analysis.

## Jupyter notebooks

We have prepared our analysis into jupyter notebooks (.ipynb files). You can either view
them on github or execute the evaluation if you install [jupyter](http://jupyter.org/). We recommend that you first
view the [Introduction.ipynb](Introduction.ipynb) file for more details.

## Installation

We recommend you install all of the depdencies through conda. Please see the miniconda [installation page](https://conda.io/miniconda.html).

Next, install the dependencies needed to run these notebooks by the following commands:

```bash
$ conda env create -f environment.yml  # create environment for CHASMplus
$ source activate UPS_jupyter  # activate environment for CHASMplus jupyter analysis
```

Remember to always activate the UPS_jupyter environment in conda!

## Runing jupyter notebooks

Using the terminal, change directories to where you download this repository. Then start jupyter lab:

```bash
$ jupyter lab
```

## Data

The notebooks use data and results available [here]().
Place the data in the top-level directory of this repository.

## Citation

We will update the citation information upon publication. 
