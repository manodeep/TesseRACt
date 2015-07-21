#########
Tutorials
#########

Running Tests
=============

Included with |tesseract| are several example halo snapshots that can be used 
to run a variety of tests. The primary function that should be used for running 
tests is :func:`tesseract.tests.run_test`. This function will create the necessary 
files, perform the tessellation, and calculate the NFW parameters using different 
techniques. For example, to test traditional fitting against tessellation based 
concentration measurement on a c = 10 halo::

    >>> tesseract.tests.run_test(c=10,nfwmeth=['fitting','voronoi'])

You can also decrease the resolution by some factor::

    >>> tesseract.tests.run_test(c=10,nfwmeth=['fitting','voronoi'],decimate=10)

or change the shape of the halo::

    >>> tesseract.tests.run_test(c=10,nfwmeth=['fitting','voronoi'],squishz=0.3)


Running Test Series
===================

In addition to individual tests, you can also run a series of tests that vary 
one parameter at a time using :func:`tesseract.tests.run_series`. For example, 
to test how tessellation based measurement of concentration is affected by how
prolate a c=50 halo is, run::

    >>> tesseract.tests.run_series('prolate',c=50,nfwmeth='voronoi')

The specific 'prolateness' values tested can be specified by either a list of values::

    >>> tesseract.tests.run_series('prolate',vlist=[0.3,0.4,0.5],c=50,nfwmeth='voronoi')

or a range and number.::

    >>> tesseract.tests.run_series('prolate',vlim=(0.3,0.5),Nv=3,c=50,nfwmeth='voronoi')


Running |tesseract| On Your Data
================================

If your data is in one of the formats currently supported by |tesseract|, 
measuring the concentration using tessellation is relatively simple. 


Creating a Parameter File
-------------------------

First you will need to create a parameter file either by hand or using 
:func:`tesseract.voro.make_param`.::

    >>> parfile = '\full\path\to\parameter\file'
    >>> tesseract.voro.make_param(parfile,**param)

where ``param`` contains all of the required keywords listed under 
:attr:`tesseract.voro._paramlist`. If you do not wish to transform or 
downsample your particles and the data does not have periodic boundary 
conditions, you can autofill many of the parameters using 
values from the provided example parameter file :attr:`tesseract.voro._example_parfile`.::

    >>> tesseract.voro.make_param(parfile,basefile=tesseract.voro._example_parfile,**param)

In this case, the only parameters you will need to provide in ``param`` are 

    * **FilePrefix**: Prefix to add to file names.
    * **PositionFile**: Full path to file containing the particle data.
    * **PositionFileFormat**: An integer specifying what format the **PositionFile** 
    is in (See :func:`tesseract.io.read_snapshot` for a list of valid formats).
    * **OutputDir**: Full path to directory where files should be saved.


Running the Tessellation
------------------------

To run the actual tessellation, use :func:`tesseract.voro.run` with the
parameter file you created in the previous step.::

    >>> tesseract.voro.run(parfile)

If you would like the output to be piped into a file rather than displayed,
provide an output file::

    >>> tesseract.voro.run(parfile,outfile='/path/to/output/file')

If the error code returned is not 0, this means that there was a problem
with the tessellation. The most common reason for a failed tessellation is 
a problem with the input snapshot. 

.. todo:: Add debug flag for voronoi?


Calculating the NFW Parameters
------------------------------

Once the tessellation is successful, NFW parameters can be computed from the 
volumes using :func:`tesseract.voro.get_nfw`::

    >>> tesseract.voro.get_nfw(parfile,method='voronoi')

Similarly, for using traditional fitting::

    >>> tesseract.voro.get_nfw(parfile,method='fitting')

or both::

    >>> tesseract.voro.get_nfw(parfile,method=['voronoi','fitting'])


Creating Snapshots
==================

If your data is in a format that is not currently supported by |tesseract|, 
fear not! The :mod:`tesseract.io` module provides utilities for writing your 
data to a format that is supported. Simply load the masses and positions for 
your particles as numpy arrays and call :func:`tesseract.io.write_snapshot`::

     >>> print type(mass), mass.shape
     (<type 'numpy.ndarray'>, (N,))
     >>> print type(pos), pos.shape
     (<type 'numpy.ndarray'>, (N,3))
     >>> filename = '\path\to\new\snapshot.dat' 
     >>> format = 0 # This is a binary file containing just masses and positions
     >>> tesseract.io.write_snapshot(filename,mass,pos,format=format)

Then just set 'PositionFile' and 'PositionFileFormat' to the new file name and 
format in your parameter file from above.