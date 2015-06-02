#!/usr/bin/python
# TODO:
# - example.param
"""
tesseract
=========

A package for measuring the concentration of halos from Nbody simulations 
non-parametrically using Voronoi tessellation.

Subpackages
-----------
voro 
  Routines for running and manipulating data returned by the Voronoi
  tesselation routine vorovol.
nfw 
  Routines relating to fitting and determining properties of NFW profiles.
io 
  Routines for data input and output.
util
  Misc. utility routines
examples
  Routines for running and plotting different tests for the provided test halos.

"""

# Basic dependencies
import ConfigParser
import os


# Config file
config_file_def = os.path.join(os.path.dirname(__file__),"default_config.ini")
config_file_usr = os.path.join(os.path.dirname(__file__),"config.ini")
config_parser = ConfigParser.ConfigParser()
config_parser.optionxform = str
config_parser.read(config_file_def)
config_parser.read(config_file_usr)
#config_parser.read(os.path.expanduser("~/.tessrc"))
#config_parser.read("config.ini")

# Fill in missing things that are necessary
#if not config_parser.has_option('general','qhulldir'):
#    qhulldir = None
#    while not isinstance(qhulldir,str) and not os.path.isdir(qhulldir):
#        config['qhulldir'] = input('Enter location of qhull installation:')
#    
# 

# General options
config = {}
if config_parser.has_option('general','outputdir'):
    config['outputdir'] = os.path.expanduser(config_parser.get('general','outputdir').strip())
else:
    config['outputdir'] = os.getcwd()

# Subpackages
import voro
import nfw
import io
import util
import examples


__all__ = ['voro','nfw','io','util','examples']


