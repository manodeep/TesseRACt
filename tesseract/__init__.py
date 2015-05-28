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
config_parser = ConfigParser.ConfigParser()
config_parser.optionxform = str
config_parser.read(os.path.join(os.path.dirname(__file__),"default_config.ini"))
config_parser.read(os.path.join(os.path.dirname(__file__),"config.ini"))
config_parser.read(os.path.expanduser("~/.tessrc"))
config_parser.read("config.ini")

# General options
config = {}
if config_parser.has_option('general','outputdir'):
    config['outputdir'] = config_parser.get('general','outputdir').strip()
else:
    config['outputdir'] = os.getcwd()

# Subpackages
import voro
import nfw
import io
import util
import examples


__all__ = ['voro','nfw','io','util','examples']


