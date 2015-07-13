--------------------------------------------------------------------------------

INTRODUCTION
============

The |tesseract| package is designed to compute concentrations of simulated dark
matter halos from volume info for particles generated using Voronoi tesselation.
This technique is advantageous as it is non-parametric, does not assume 
spherical symmetry, and allows for the presence of substructure. For a more
complete description of this technique including a comparison to other 
techniques for calculating concentration, please see the accompanying paper 
`Lang et al. (2015)`_.

This package includes:

* **vorovol**: C program for computing the Voronoi diagram of particle data in 
    a number of formats including Gadget-2, Gasoline, binary, and ASCII as well
    as BGC2 halo catalogues.
* routines for compiling, running, and parsing **vorovol** output
* routines for computing concentrations using particles volumes, traditional 
    fitting to an NFW profile, and non-parametric techniques that assume 
    spherical symmetry
* routines and test halos for running many of the performance tests presented in 
    `Lang et al. (2015)`_

Some useful links associated with |tesseract|:
* `PyPI`_ - The most recent stable release.
* `Docs`_ - Tutorials and descriptions of the package modules and functions.
* `Lang et al. (2015)`_ - The scientific paper associated with this package.

If you would like more information about |tesseract|, please contact `Meagan Lang`_.

.. _PyPI: https://pypi.python.org/pypi/tesseract                                                                                                                                      
.. _Lang et al. (2015): http://arxiv.org/abs/1504.04307                                                                                                                               
.. _Meagan Lang: cfh5058@gmail.com                                                                                                                                                    
.. _qhull: http://www.qhull.org/                                                                                                                                                      
.. |tesseract| replace:: TesseRACt  
.. _Docs: https://readthedocs.org/projects/pytesseract/