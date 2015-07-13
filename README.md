--------------------------------------------------------------------------------

INTRODUCTION
============

The tesseract package is designed to compute concentrations of simulated dark
matter halos from volume info for particles generated using Voronoi tesselation.
This technique is advantageous as it is non-parametric, does not assume 
spherical symmetry, and allows for the presence of substructure. For a more
complete description of this technique including a comparison to other 
techniques for calculating concentration, please see the accompanying paper 
Lang et al. (2015).

This package includes:

* **vorovol**: C program for computing the Voronoi diagram of particle data in 
    a number of formats including Gadget-2, Gasoline, binary, and ASCII as well
    as BGC2 halo catalogues.
* routines for compiling, running, and parsing **vorovol** output
* routines for computing concentrations using particles volumes, traditional 
    fitting to an NFW profile, and non-parametric techniques that assume 
    spherical symmetry
* routines and test halos for running many of the performance tests presented in 
    Lang et al. (2015)
    
Full documentation for TesseRACt is hosted on [readthedocs.org](https://readthedocs.org/projects/pytesseract/).


#TODO
#====
#Things that are being worked on:
#
#* create example parameter file
#* create straight forward tests
#* Add runner test script for install
#* add outputfile parameter to vorovol