import os,pickle
import numpy as np
from . import util
from . import config,config_parser,_config_file_usr
# TODO: 
# - fix units in snapshots and volume files so that Mscl and Rscl are not
#   needed in get_nfw
# - suppress vorovol output from qhull

_installdir = os.path.dirname(os.path.realpath(__file__))
_installdir_vorvol = os.path.join(_installdir,'vorovol')

# Vorovol files
_makefile_vorovol = os.path.join(_installdir_vorvol,'Makefile')
_execfile_vorovol = os.path.join(_installdir_vorvol,'vorovol')

# Qhull installation
if config_parser.has_option('qhull','install-dir'):
    _qhulldir = config_parser.get('qhull','install-dir').strip()
else:
    # Ask user for qhull directory
    qstr = 'Please enter the path to an existing directory where qhull should be installed: '
    _qhulldir = os.path.expanduser(raw_input(qstr).strip())
    while not os.path.isdir(_qhulldir):
        print 'That is not a valid directory.'
        _qhulldir = os.path.expanduser(raw_input(qstr).strip())
    # Add option to config file
    if not config_parser.has_section('qhull'):
        config_parser.add_section('qhull')
    config_parser.set('qhull','install-dir',_qhulldir)
    with open(_config_file_usr,'w') as fp:
        config_parser.write(fp)
    # Unpack the tar file
    _qhulltar = os.path.join(_installdir,'qhull2002.1.tar')
    if os.path.isdir(os.path.join(_qhulldir,'qhull2002.1')):
        print 'There is already a qhull installation there. No new installation necessary.'
    else:
        print 'Unpacking Qhull in {}...'.format(_qhulldir)
        os.system('tar -C {} -xf {}'.format(_qhulldir,_qhulltar))

# Qhull files
_installdir_qhull = os.path.join(_qhulldir,'qhull2002.1','src')
_makefile_qhull = os.path.join(_installdir_qhull,'Makefile')
_execfile_qhull = os.path.join(_installdir_qhull,'qhull_a.h')
os.environ['QHULLSRCDIR'] = _installdir_qhull

# Parameter file options
_paramlist = ['FilePrefix','FileSuffix','NumDivide','PeriodicBoundariesOn',
              'Border','BoxSize','PositionFile','PositionFileFormat','OutputDir',
              'OutputAdjacenciesOn','MaxNumSnapshot',
              'DecimateInputBy','SquishY','SquishZ',
              'GadgetParticleType','BgTreebiNskip','Bgc2HaloId']
_paramopt = ['NumDivide','OutputAdjacenciesOn','MaxNumSnapshot',
             'DecimateInputBy','SquishY','SquishZ',
             'GadgetParticleType','BgTreebiNskip','Bgc2HaloId']

# ------------------------------------------------------------------------------
# METHODS FOR INTERFACING WITH THE VOROVOL ROUTINE
def run(parfile0,exefile=None,outfile=None,overwrite=False,verbose=True,
        recompile=False):
    """
    Runs vorovol on the designated parameter file, compiling first if necessary.
    A non-zero return code refers to an error generated while running vorovol.
        parfile  : path to vorovol parameter file or a dictionary of parameters
                   that should be written to a file. Note: one of the keys can
                   be 'parfile' and contain the parameter file path. If it is
                   not, a generic parameter file named 'voro.param' is created.
        exefile  : path to vorovol executable file
        outfile  : path to file where output from vorovol should be piped. If 
                   not provided, output will be printed to the screen.
        overwrite: if true, the executable is run even if output for this
                   run already exists. (default = False)
        verbose  : if true, the run prints out information. (default = True)
        recompile: if true, the executable is compiled even if it exists.
                   (default = True)
    """
    # Set defaults and allow for parfile or parameters
    if exefile is None: exefile = _execfile_vorovol
    if isinstance(parfile0,dict):
        param = parfile0
        parfile = param.get('parfile','voro.param')
        if not os.path.isfile(parfile):
            param = make_param(parfile,**param)
    else:
        parfile = parfile0
        if not os.path.isfile(parfile):
            raise Exception('Parameter file does not exist: {}'.format(parfile))
        else:
            param = read_param(parfile)
    # Check for output
    volfile = namefile('vols',param)
    if os.path.isfile(volfile) and not overwrite:
        if verbose:
            print 'voro.py @ 63: Voronoi output already exists and overwrite not set.'
            print '    '+volfile
        return 0
    # Compile executable if it does not exists
    if not os.path.isfile(exefile) or recompile:
        make_vorovol(os.path.join(os.path.dirname(exefile),'Makefile'))
    # Create output directories
    outputdir = param['OutputDir']
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    for d in ['vols','adjs','part']:
        if not os.path.isdir(os.path.join(outputdir,d)):
            os.mkdir(os.path.join(outputdir,d))
    # Create command
    cmd = './{} {}'.format(os.path.basename(exefile),parfile)
    if outfile is not None:
        cmd+=' > {}'.format(outfile)
    # Run
    curdir = os.getcwd()
    os.chdir(os.path.dirname(exefile))
    if verbose:
        print cmd
    code = os.system(cmd)
    os.chdir(curdir)
    # Return code
    return code


def make(makefile,exefile=None):
    """
    Compiles an application using a Makefile.
        makefile: path to Makefile
    """
    # Set default makefile & check that it exists
    makefile = os.path.expanduser(makefile)
    if not os.path.isfile(makefile):
        raise Exception('Makefile does not exists: {}'.format(makefile))
    # Record current working directory then move to executable directory
    curdir = os.getcwd()
    os.chdir(os.path.dirname(makefile))
    # Clean up if the exeutable exists
    os.system('make clean')
    # Make & check that it worked
    os.system('make')
    if isinstance(exefile,str) and not os.path.isfile(exefile):
        raise Exception('Make failed to created executable: {}'.format(exefile))
    # Move back to original directory
    os.chdir(curdir)
    return

def make_vorovol(makefile=_makefile_vorovol,makefile_qhull=_makefile_qhull):
    """
    Compiles vorovol executable using a Makefile.
    """
    # Check for required qhull library
    execfile_qhull = os.path.join(os.path.dirname(makefile_qhull),'qhull_a.h')
    if not os.path.isfile(execfile_qhull):
        make_qhull(makefile_qhull)
    make(makefile,exefile=os.path.join(os.path.dirname(makefile),'vorovol'))
    return

def make_qhull(makefile=_makefile_qhull):
    """
    Compiles qhull libraries using a Makefile.
    """
    make(makefile,exefile=os.path.join(os.path.dirname(makefile),'qhull_a.h'))
    return

# ------------------------------------------------------------------------------
# METHODS FOR HANDLING VOLUME FILES
def namefile(name,param):
    """
    Returns filename of voronoi output file
        name : string specifying the type of file (part,vols,adjs)
               Note: Adjacencies will only be output if OutputAdjacenciesOn 
               is 1 in the parameter file.
        param: vorovol parameters
    """
    # Parse name
    listnames = ['part','vols']
    if param['OutputAdjacenciesOn']:
        listnames.append('adjs')
    else:
        if name=='adjs':
            raise Exception('The output from this run does/will not contain adjacencies. \n'+
                            'Set OutputAdjacenciesOn to 1 in the parameter file.')
    if name not in listnames:
        raise Exception('Valid output file names include: {}'.format(listnames))
    # Base filename and directory
    fdir = os.path.join(param['OutputDir'],name)
    fnam = param['FilePrefix']
    # Add strings
    if param['DecimateInputBy']>1: 
        fnam+='_dec{:d}'.format(param['DecimateInputBy'])
    if param['PositionFileFormat']==1 and param['GadgetParticleType']>=0: 
        fnam+='_{:d}'.format(param['GadgetParticleType'])
    if param['PositionFileFormat']==3 and param['Bgc2HaloId']>=0:
        fnam+='_{:d}'.format(param['Bgc2HaloId'])
    if param['SquishY']>0 and param['SquishY']<1:
        fnam+='_sqY{:4.2f}'.format(param['SquishY']).replace('.','p')
    if param['SquishZ']>0 and param['SquishZ']<1:
        fnam+='_sqZ{:4.2f}'.format(param['SquishZ']).replace('.','p')
    # Extensions
    if name == 'part':
        ext = '.*' # These files have 00.00.00 format, numbers depend on split
    else:
        ext = ''
    return os.path.join(fdir,'{}{}.{}{}'.format(fnam,param['FileSuffix'],name,ext))

def read_snapshot(param,return_npart=False,**kwargs):
    """
    Read position file. Downsample and performing volume preserving 
    transformation based on parameters. (This replicates what vorovol does when
    it loads snapshots.)
        param       : Dictionary of vorovol parameters or path to parameter 
                      file containing them.
        return_npart: Set to True if only the number of particles in the 
                      snapshot should be returned. (default = False)
    Additional keywords are passed to the appropriate method for reading the
    PositionFileFormat specified in the parameters.
    """
    from . import io
    # Read parameters if its a string
    if isinstance(param,dict): pass
    elif isinstance(param,str):
        param = read_param(param)
    else:
        raise Exception('Invalid parameter type: {}.'.format(type(param))+
                        'Must be a dictionary or path to a parameter file.')
    # Get necessary parameters
    filename = param['PositionFile']
    format = param['PositionFileFormat']
    # Get format specify keywords for the read function
    if format in io._snapshot_formats:
        kwargs.update(**io._snapshot_formats[format].parse_voropar(param))
    else:
        raise Exception('There is not a snapshot type associated with the '+
                        'format code {} in io.py'.format(format))
    # Load masses and positions
    kwargs.update(format=format,return_npart=return_npart)
    out = io.read_snapshot(filename,**kwargs)
    # Parse output
    if return_npart:
        Nold = out
    else:
        mass,pos = out
        Nold = len(mass)
    # Determine number of particles after down-sampling
    if param['DecimateInputBy']>1:
        Nnew = Nold/param['DecimateInputBy']
    else:
        Nnew = Nold
    # Return npart
    if return_npart:
        return Nnew
    # Squeeze
    pos = util.squeeze_constvol(pos,yfact=param['SquishY'],
                                zfact=param['SquishZ'])
    # Decimate particle number and increase mass of particles to compensate
    if param['DecimateInputBy']>1:
        idxdec = slice(0,Nold - (Nold % param['DecimateInputBy']),
                       param['DecimateInputBy'])
        pos = pos[idxdec,:]
        mass = mass[idxdec]
        if len(mass)!=Nnew:
            raise Exception('Error downsampling. Should have {}, '.format(Nnew)+
                            'but we have {}.'.format(len(mass)))
        mass*= (float(Nold)/float(Nnew))
    # Return
    return mass,pos

def read_volume(filename):
    """Reads vorovol volume file"""
    import struct
    # Read in volumes
    fd = open(filename,"rb")
    npart = struct.unpack('I',fd.read(4))[0]
    vols = np.fromfile(fd,dtype=np.float32,count=npart)
    # Close and return
    fd.close()
    return vols

def write_volume(filename,vols,overwrite=False):
    """Writes vorovol volume file"""
    import struct
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        return
    # Write volumes
    fd = open(filename,"wb")
    fd.write(struct.pack('I',len(vols)))
    vols.tofile(fd)
    # Close and return
    fd.close()
    return
    

# ------------------------------------------------------------------------------
# METHODS FOR INTERFACING WITH PARAMETER FILES
def make_param(filename,basefile=None,overwrite=False,**kwargs):
    """
    Create vorovol parameter file
        filename : path to file where parameters should be saved
        basefile : path to parameter file that new file should be created from
                   (default = None, file is created from scratch and all 
                   required fields must be provided)
        overwrite: if True, existing file is replaced (default = False)
        Additional keywords are parameter fields including:
        FilePrefix          : Prefix to add to file names
        FileSuffix          : Suffix to add to end of file names (this can be a
                              C style format code if position file refers to 
                              multiple snapshots)
        NumDivide           : Number of times snapshot is divided
        PeriodicBoundariesOn: Set to 1 if snapshot is periodic
        Border              : Border added to edge of boxes
        BoxSize             : Size of periodic box
        PositionFile        : Path to file containing positions (this can 
                              contain a C style format code in order to refer
                              to multiple snapshot with the same naming scheme)
        PositionFileFormat  : Integer specifying position file format
         -1: Pre-existing vorovol parts file
          0: f77 unformatted binary snapshot
          1: Gadget snapshot
          2: Buildgal TREEBI files
          3: BGC2 halo catalogue
        GadgetParticleType  : Integer gadget particle type
         -1: all particles
          0: gas particles
          1: dark matter particles
          2: disk particles
          3: bulge particles
          4: star particles
          5: boundary particles, but why?
        BgTreebiNskip       : Number of particles to skip in TREEBI snapshot
        Bgc2HaloId          : ID number of halo in Bgc2 catalogue
        OutputDir           : Path to directory where output should be saved
        OutputAdjacenciesOn : Set to 1 to also output adjacencies info
        DecimateInputBy     : >1 factor to decimate particle number by
        MaxNumSnapshot      : Maximum number of snapshot to do (only invoked if
                              FileSuffix contains a format code)
        SquishY             : factor to squish Y coords (preserves volume)
        SquishZ             : factor to squish Z coords (preserves volume)
    """
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Read base file if provided
    if basefile is not None:
        param = read_param(basefile)
    else:
        param = {}
    # Fill in values from keywords
    for k in _paramlist:
        if k in kwargs: 
            param[k] = kwargs[k]
    # Create directory
    pdir = os.path.dirname(filename)
    if not os.path.isdir(pdir):
        os.mkdir(pdir)
    # Write to file
    write_param(filename,param,overwrite=overwrite)
    return param

def write_param(filename,param,overwrite=False):
    """
    Writes vorovol parameter file
        filename : path to file where parameters should be written
        param    : dictionary of parameters to write
        overwrite: if True, existing file is replaced (default = False)
    """
    import copy
    width=29
    optpar = copy.deepcopy(_paramopt)
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Add optional parameters
    if param.get('PeriodicBoundariesOn',1)==0:
        optpar+=['BoxSize','Border']
    # Check for missing parameters
    missing = []
    for k in _paramlist:
        if k not in param:
            if k in optpar:
                param[k] = -1
            else:
                missing.append(k)
    if len(missing)>0:
        raise Exception('There are missing parameter fields: {}'.format(missing))
    # Open file
    fd=open(filename,'w')
    # Loop over parameters
    for k in _paramlist:
        v=param[k]
        fd.write(k.ljust(width)+util.val2str(v)+'\n')
    # Close file & return
    fd.close()
    return

def read_param(filename):
    """
    Read vorovol parameter file. Parameter are returned in a dictionary.
        filename: file parameters should be read from
    """
    cchar = '%'
    # Open and initialize
    fd = open(filename,'r')
    param={}
    keylist=[]
    # Loop over lines
    for line in fd:
        # Skip comments
        if line.startswith(cchar):
            pass
        else:
            # Remove comments and split key/value
            line=line.split(cchar)[0]
            vars=line.split()
            if len(vars)==2:
                key=vars[0]
                val=util.str2val(vars[1])
            else:
                key=vars[0]
                val=util.str2val('')
            param[key]=val
            keylist.append(key)
    # Close
    fd.close()
    # Check for missing parameters
    missing = []
    for k in _paramlist:
        if k not in param:
            missing.append(k)
    if len(missing)>0:
        raise Exception('There are missing parameter fields: {}'.format(missing))
    # Return
    return param


# ------------------------------------------------------------------------------
# METHODS TO COMPUTE NFW
def get_nfw(param,method='voronoi',vorometh='rhalf',nfwfile=None,
            ownfw=False,plotflag=False,plotfile=None,residuals=True,delta=None,
            Mscl=1.,Rscl=1.,**kwargs):
    """
    Returns NFW parameters found using the specified method(s). The output
    is a dictionary of NFW parameters. If more than one method is specified, 
    the returned dictionary has the listed methods as keys with the 
    corresponding NFW parameter dictionaries as values.
        param   :  dictionary of vorovol parameters or path to a parameter file
        method  :  method or list of methods that should be used to calculate NFW
                   parameters. See nfw.calc_nfw for more detailed info on 
                   techniques (default = 'voronoi')
          'voronoi'   : use unweighted voronoi volumes to compute radii (See 
                        util.vol2rad) and then use technique specified by
                        vorometh.
          'voronoi_wX': same as 'voronoi', but weight voronoi volumes by
                        the geometric volumes computed from the actual radii.
                        The weight X specifies how much the geometric volumes
                        should be weighted. 
          'fit'       : leastsq fit to NFW enclosed mass profile
          'rhalf'     : non-parametric half-mass
          'vpeak'     : non-parametric peak velocity
          Any valid inputs for the method argument of nfw.calc_nfw are also 
          valid here.
       vorometh:  method that should be used to compute concentration using
                  voronoi based radii. This is only used for the 'voronoi' and
                  'voronoi_wX' methods. Valid values include 'fit', 'rhalf',
                  and 'vpeak' (See above).
       nfwfile  : path to file where NFW data should be saved. If None, the data
                  is just returned.
       ownfw    : if True, existing nfwfile is overwritten. (default = False)
       plotflag : if True, the data and profile are plotted (default = False)
       plotfile : path to file where plot should be saved. If None, the plot
                  is displayed instead. (default = None)
       residuals: if True, residuals are also plotted
       delta    : virial overdensity factor (default = 200 for spherical 
                  techniques, 243 for voronoi)
       Additional keywords are passed to each calc_nfw method that is called.
    """
    from . import nfw
    colors = ['b','r','g','m','c']
    limdef = (-0.01,0.01)
    # Create list of methods
    if isinstance(method,list):
        methlist = method
    else:
        methlist = [method]
    # Read NFW file if it exists
    if isinstance(nfwfile,str) and os.path.isfile(nfwfile) and not ownfw:
        fd = open(nfwfile,'r')
        out = pickle.load(fd)
        fd.close()
        miss = False
        for m in methlist:
            if m not in out:
                miss = True
        if not miss:
            if isinstance(method,list):
                return out
            else:
                return out[method]
    else:
        out = {}
    # Read parameters
    if isinstance(param,str):
        param = read_param(param)
    elif isinstance(param,dict):
        pass
    else:
        raise Exception('param must be dictionary or path to parameter file '+
                        'not {}'.format(type(param)))
    # Read position file
    mass,pos = read_snapshot(param)
    mass*=Mscl
    pos*=Rscl
    rad = np.sqrt(pos[:,0]**2. + pos[:,1]**2. + pos[:,2]**2.)
    radsort = np.argsort(rad)
    # Setup plot
    if plotflag:
        if len(colors)<len(methlist):
            raise Exception('Only {} colors for {} methods.'.format(len(colors),len(methlist)))
        import matplotlib.pyplot as plt
        plt.clf()
        fig = plt.figure(figsize=(10,5))
        if residuals:
            axs = fig.add_axes((.1,.3,.6,.6))
            plt.setp( axs.get_xticklabels(), visible=False)
            axs_res = fig.add_axes((.1,.1,.6,.2))
            axs_res.set_xlim(limdef)
            axs_res.set_ylim(limdef)
        else:
            axs = plt.subplot(1,1,1)
            axs_res = None
        axs.set_xlim(limdef)
        axs.set_ylim(limdef)
    else:
        axs = None
        axs_res = None
    # Loop over method
    vol = None
    for i,m in enumerate(methlist):
        if m in out:
            continue
        # Voronoi based method
        if m.startswith('voronoi'):
            if delta is None:
                idelta = 200.#243.
            else:
                idelta = delta
            # Read volumes and remove infinite values (negatives)
            if vol is None:
                vol = read_volume(namefile('vols',param))
                vol*=(Rscl**3.)
                idxfin = (vol>=0)
            # Get weights
            if m.startswith('voronoi_w'):
                weight = float(m.split('voronoi_w')[-1])
                weightby = rad[idxfin]
            else:
                weight = False
                weightby = False
            # Compute radii and sorting
            volrad,idxsort = util.vol2rad(vol[idxfin],outsort=True,
                                          weightby=weightby,weight=weight)
            # Virial things
            #rho = mass[idxfin][idxsort]/vol[idxfin][idxsort]
            # rvir,mvir = nfw.calc_virial(volrad,rho=rho,delta=idelta,
            #                             rhoc=kwargs.get('rhoc',None))
            #print 'voro.py @ 424: rvir,mvir = ',rvir,mvir
            # Get nfw params
            out[m] = nfw.calc_nfw(volrad,m=mass[idxfin][idxsort],
                                  method=vorometh,issorted=True,
                                  plotflag=plotflag,residuals=residuals,
                                  axs=axs,axs_res=axs_res,label=m.title(),
                                  color=colors[i],delta=idelta,**kwargs)
            # Scale based on voronoi radii calibration
            #fit_vor2rad = [1.20765393, 0.96956452]
#            fit_vor2rad = [ 1./1.24032053,  1./0.95997943]
            fit_vor2rad = [ 0.84497027,  1.02418653]
            c = out[m]['c']
            out[m]['c'] = fit_vor2rad[0]*(c**fit_vor2rad[1])
        # Radial based methods
        else:
            if delta is None:
                idelta = 200.
            else:
                idelta = delta
            out[m] = nfw.calc_nfw(rad[radsort],m=mass[radsort],method=m,
                                  issorted=True,plotflag=plotflag,residuals=residuals,
                                  axs=axs,axs_res=axs_res,label=m.title(),
                                  color=colors[i],delta=idelta,**kwargs)
        out[m]['N'] = len(mass)
    # Save data
    if isinstance(nfwfile,str):
        fd = open(nfwfile,'w')
        pickle.dump(out,fd)
        fd.close()
    # Save/show plot
    if plotflag:
        # Legend
        leg = axs.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # Remove bottom ticks
        if residuals:
            fig.canvas.draw()
            tics = [item.get_text() for item in axs.get_yticklabels()]
            tics[0] = ''
            axs.set_yticklabels(tics)
        # Save/show
        if plotfile is None:
            plt.show()
        else:
            plt.savefig(plotfile,bbox_inches='tight',bbox_extra_artists=[leg])
            print '    '+plotfile
    # Return output
    if not isinstance(method,list):
        return out[method]
    else:
        return out
