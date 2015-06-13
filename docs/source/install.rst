############
Installation
############

|tesseract| can be installed from either `PyPI`_ or from the source distribution.

Installing from PyPI
====================

The easiest way to install |tesseract| is using ``pip``. If you have administrative privleges on the target machine, this is done using::

   $ pip install tesseract

If you do not have admin privleges, simply install it locally using::

   $ pip install tesseract --user

The |tesseract| package can then be updated to the most recent stable release using::

   $ pip install tesseract --upgrade

Installing from the Source Distribution
=======================================

The most recent |tesseract| source distribution can be obtained by either downloading or cloning the repository from Bitbucket. Using `Mercurial <https://mercurial.selenic.com/>`_ this is done by issuing the following command::

   $ hg clone https://[username]@bitbucket.org/[username]/tesseract

where ``[username]`` should be replaced with your Bitbucket username. The Bitbucket repository is currently private. If you would like access to this repository, please contact `Meagan Lang`_. 

.. todo:: Add explanation on cloning the repository using git, svn, etc.

Once you have the |tesseract| source distribution, move into the distribution directory::

   $ cd tesseract

and use the standard `Distutils <https://docs.python.org/2/distutils/>`_ command to build and install the distribution::

   $ python setup.py install

If you do not have administrative privleges, this can be done using::

   $ python setup.py install --user 

Testing the Install
===================

To test that everything was installed propertly. From the python prompt, import |tesseract|::

   >>> import tesseract

and try to access the documentation::

   >>> help(tesseract)

Additional tests can be found in :doc:`tests`.

The First Import
================

The first time you  run	``import tesseract``, a few things will happen. First, a user config file ``.tessrc`` will be created in your home directory. This file is used to control different aspects of |tesseract| which are explained in :doc:`config`. Second, you will be prompted to enter a directory in which `qhull`_ will be installed. This directory will be added to the user configuration file, which can be changed at any time if you move `qhull`_.


     



