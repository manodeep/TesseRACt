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

* vorovol: a c program for computing the Voronoi tessellation for particle data in a number of formats including Gadget-2 snapshots, ascii, and binary files.
* routines for compiling, running, and parsing vorovol output
* routines for computing concentrations using particles volumes, traditional fitting to an NFW profile, and non-parametric techniques that assume 
* routines and test halos for running many of the performance tests presented in Lang et al. (2015)

--------------------------------------------------------------------------------

INSTALLATION
============

Only installation through cloning the source code is currently supported.

The easiest way to install TesseRACt from the source code is to clone the repository using mercurial. 

1. **Install mercurial if you don't have it.** This is easy to do via `pip install hg` or `pip install hg --user` if you do not have admin privileges. 
2. **Clone the TesseRACt repository.** Move to the directory where you want to install TesseRACt. Then issue the command `hg clone [https link]` where `[https link]` should be the bitbucket link.

--------------------------------------------------------------------------------

TODO
====
Things that are being worked on:

* liscense
* create example parameter file
* customize Makefile on install
* ask user if qhull should be downloaded & update parameter file
* create straight forward tests
* create config file containing output dir
* add outputfile parameter to vorovol