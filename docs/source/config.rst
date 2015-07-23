###############
The Config File
###############

The |tesseract| user config file `.tessrc` is created in your home directory when |tesseract| is first imported. It is initialized from the default configuration file `default_config.ini` and can be edited at any time to change different |tesseract| aspects. However, the default configuration file should NOT be edited directly in case new functionality is added. Any options not provided in `.tessrc` are initialized to the values in `default_config.ini`. 

For additional details on configuration file syntax, please see the documentation for the `ConfigParser <https://docs.python.org/2/library/configparser.html>`_ package.

General Options
===============

**outputdir**
   Path to directory where |tesseract| output should be saved. If not provided, output will be saved in the current working directory.


NFW Options
===========

Options for controlling how NFW profile parameters are computed.

**default-rhoc**
   The default value used for the critical density of the universe in units of `Msol kpc**-3`. (`1.1845e2 Msol kpc**-3` by default)

**default-delta**
   The default value used to define the virial over-density. This is used to calculate the virial radius and mass of a halo. (200 by default)

Test Options
============

Options for controlling the examples provided.

**snapshot-format**
   Code specifying the format of the test halo snapshots. This should not be edited unless you convert the test halo snapshots into another snapshot format. (0 by default)
**halodir**
   Directory containing test halo snapshots. If not provided, it is assumed to be the directory within the |tesseract| distribution.
**copydir** 
   Directory containing copies of the test halo snapshots initialized with different random number seeds. As these snapshots are not provided with the public distribution, this option should not be used.
**avail-conc**
   List of concentrations of the test halo snapshots. The current version of |tesseract| includes test halos with concentrations of 5, 10, 25, and 50.
**default-series** 
   String specifying test series that should be run by default when `examples.run_series` is called. ('conc' by default)
**avail-series** 
   List of available test series. |tesseract| currently supports the following series which are described in the Examples section below: conc, oblate, prolate, triax, npart, substr_mass, substr_rsep, substr_conc, and substr_rho0.
**nfw-methods** 
   List of methods for calculating NFW parameters that are used for each test by default. (voronoi, fit, rhalf, and vpeak by default)
**default-conc** 
   Default concentration used for tests (10 by default)
**default-subm** 
   Default mass of subhalo used for substructure tests in terms of the parent halo's virial mass (0.1 by default)
**default-subr**
   Default radius that subhalo is placed at for substructure tests in terms of the parent halo's virial radius (0.5 by default)
**default-subc**
   Default concentration of subhalo used for substructure tests. (50 by default)
**default-subrho**
   Default density of subhalo used for substructure tests in terms of the parent halo's central concentration (0.5 by default)
**default-ellip**
   Default ellipticity of halo for test which vary triaxiality. (0.5 by default)

qhull Options
=============

Options for controlling `qhull`_.

**install-dir** 
   The directory under which qhull should be installed. This is initialized the first time that the |tesseract| package is imported. If this does not point to a valid qhull installation, |tesseract| will be unable to perform the Voronoi tessellation.



