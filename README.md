# How to use the sterile-neutrinos-constraints repository
This repository complements the paper "Constraints on sterile neutrino models from strong gravitational lensing, Milky Way
satellites, and Lyman - α forest", by Zelko et al 2022 ,Phys. Rev. Lett. **129**, 191301
https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.129.191301

### Library Requirenments


To run the python code,  you will need the following packages installed:

matplotlib, numpy, pickle, os, scipy, sys,

as well as jupyter notebook.

I made the code with python version 3.8.5, and hopefully it should work with other versions.

### System environment variables
You will have to set two environment variables in your operating system:

DARK_MATTER_PAPER_LOCATION - a path to where you would like the paper plots and tables to be saved by the scripts.
DARK_MATTER_CODE_LOCATION - a path to the directory where you downloaded the code for this project.

This [guide](https://www.twilio.com/blog/2017/01/how-to-set-environment-variables.html) explains how to set environment variables in Windows, MacOS, Linux.


### How to reproduce the results of the paper

Please see the **paper_analysis.ipynb** notebook to see how the paper results are obtained.

The Mathematica folder contains the mathematica notebook that shows how the relations in the Appendix are obtained.

The strong_lensing_data folder contains the pickled python function that returns the 2D likelihood function from [Gilman et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.6077G/abstract) used in the analysis; the python file strong_lensing_data_processing.py contains the functions to process this.

the transfer_function folder contains the transfer functions for the Higgs decay model from [Petraki & Kusenko 2008 (PK)](https://ui.adsabs.harvard.edu/abs/2008PhRvD..77f5014P/abstract) and the GUT scale scenario from [Kusenko et al 2010 (KTY)](https://ui.adsabs.harvard.edu/abs/2010PhLB..693..144K/abstract), as obtained by [Abazajian & Kusenko 2019](https://arxiv.org/abs/1907.11696). The python script transfer_functions.py has the functions to process these.

The other data used in the paper are taken directly from their respective papers.
