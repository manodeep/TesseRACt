from distutils.core import setup
setup(
    name = 'tesseract',
    packages = ['tesseract'],
    package_dir = {'tesseract': 'tesseract'},
    package_data = {'tesseract': ['README.txt','default_config.ini','halos/*.snap',
                                  'vorovol/*.c','vorovol/*.h','vorovol/Makefile']},
    version = '0.1.0',
    author = 'Meagan Lang',
    author_email = 'meagan.lang@vanderbilt.edu',
    url = 'http://vpac00.phy.vanderbilt.edu/~langmm/index.html',
    description = 'Tesselation based Recovery of Amorphous halo Concentrations',
    classifiers = ["Programming Language :: Python","Operating System :: OS Independent",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Intended Audience :: Science/Research","Natural Language :: English",
                   "Topic :: Scientific/Engineering","Topic :: Scientific/Engineering :: Astronomy",
                   "Development Status :: 3 - Alpha"], #"Programming Language :: C" 
    long_description = """\
Tesselation based Recovery of Amorphous halo Concentrations
-----------------------------------------------------------

Computes concentrations of simulated dark matter halos non-parametrically 
particle associated volumes.


"""
)
