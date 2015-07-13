============
Introduction
============

The tesseract package is designed to compute concentrations of simulated dark
matter halos from volume info for particles generated using Voronoi tesselation.
This technique is advantageous as it is non-parametric, does not assume 
spherical symmetry, and allows for the presence of substructure. For a more
complete description of this technique including a comparison to other 
techniques for calculating concentration, please see the accompanying paper 
`Lang et al. (2015)`_.

This package allows users to:

* perform Voronoi tessellation through access to `qhull`_ routines
* measure particle distribution properties like concentration using different techniques including tessellation
* replicate the performance test presented in `Lang et al. (2015)`_



