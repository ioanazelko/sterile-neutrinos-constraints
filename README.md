# How to use the sterile-neutrinos-constraints repository
This repository complements the paper "Constraints on sterile neutrino models from strong gravitational lensing, Milky Way
satellites, and Lyman - α forest", by Zelko et al 2022 (currently in the journal review stage)

### Library Requirenments


To run the python code,  you will need the following packages installed:

matplotlib, numpy, pickle, os, scipy, sys,

as well as jupyter notebook.

I made the code with python version 3.8.5, and hopefully it should work with other versions.

### System environment variables
You will have to set two environment variables in your operating system:

DARK_MATTER_PAPER_LOCATION - a path to where you would like the paper plots and tables to be saved by the scripts
DARK_MATTER_CODE_LOCATION - a path to the directory where you downloaded the code for thir project

This is a guide explaining how to set environment variables in Windows, MacOS, Linux
https://www.twilio.com/blog/2017/01/how-to-set-environment-variables.html

### How to reproduce the results of the paper

Please see the paper_analysis.pybn notebook to see how the paper results are obtained.

The Mathematica folder contains the mathematica notebook that shows how the relations in the Appendix are obtained.

The strong_lensing_data folder contains the pickled python function that returns the 2D likelihood function from Gilman et al. 202 used in the analysis; the python file strong_lensing_data_processing.py contains the functions to process this.

the transfer_function folder contains the transfer functions for the Higgs decay model from Petraki & Kusenki 2008 (PK) and the GUT scale scenario from Kusenki et al 2010 (KTY), as obtained by Abazajian & Kusenko 2019. The python script transfer_functions.py has the functions to process these.
