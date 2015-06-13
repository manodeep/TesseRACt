import numpy as np
import os

# TODO:
# - Gadget IO

# ------------------------------------------------------------------------------
# LOGGING OF SNAPSHOT FORMATS
_snapshot_formats = {}
def register_snapshot_format(code):
    """Decorator to register snapshot format classes"""
    def wrapper(f):
        if not issubclass(f,Snapshot):
            raise Exception('Only Snapshot subclasses can be registered.')
        if code in _snapshot_formats:
            raise Exception('Snapshot format code {} was already '.format(code)+
                            'assigned to {}. '.format(_snapshot_formats[code])+
                            'Cannot also be assigned to {}.'.format(f.__name__))
        f._code = code
        _snapshot_formats[code] = f()
        _snapshot_formats[code]._code = code
        _snapshot_formats[code].__name__ = f.__name__
        return f
    return wrapper
def display_snapshot_formats():
    """Prints information on the registered snapshot formats"""
    print 80*'='
    print '{:4s}  {:20s}  {}'.format('code','name','description')
    print 80*'-'
    for c in sorted(_snapshot_formats.keys()):
        f = _snapshot_formats[c]
        print '{:4d}  {:20s}  {}'.format(c,f.__name__,f.__doc__)
    print 80*'='

class Snapshot(object):
    """Base class for snapshot formats"""
    @property
    def code(self):
        if not hasattr(self,'_code'):
            raise AttributeError('This snapshot type does not have a code.')
        return self._code
    def read(self,*args,**kwargs):
        """Dummy read method to return error"""
        raise Exception('This snapshot type does not have a read method.')
    def write(self,*args,**kwargs):
        """Dummy write method to return error"""
        raise Exception('This snapshot type does not have a write method.')
    def load(self,*args,**kwargs):
        """Alias for read."""
        return self.read(*args,**kwargs)
    def parse_voropar(self,param):
        """
        Default method for parsing vorovol parameters. Returns an empty
        dictionary.
        """
        return {}

# ------------------------------------------------------------------------------
# GENERAL SNAPSHOT FILE HANDLING
def read_snapshot(filename,format=0,**kwargs):
    """
    Read masses and positions from a snapshot.
        filename: Full path to the file that should be read.
        format  : Integer specifying the type of snapshot. See output from
                  display_snapshot_formats for a description of the available 
                  codes. (default = 0)
    Additional keywords are passed to the appropriate method for reading.
    """
    if format in _snapshot_formats:
        out = _snapshot_formats[format].read(filename,**kwargs)
    else:
        raise Exception('No snapshot format registered with code '+
                        '{}.'.format(format))
    return out

def write_snapshot(filename,mass,pos,format=0,**kwargs):
    """
    Write masses and positions to a snapshot.
        filename: Full path to the file that should be written.
        format  : Integer specifying the type of snapshot. See output from
                  display_snapshot_formats for a description of the available 
                  codes. (default = 0)
    Additional keywords are passed to the appropriate method for writing.
    """
    if format in _snapshot_formats:
        out = _snapshot_formats[format].write(filename,mass,pos,**kwargs)
    else:
        raise Exception('No snapshot format registered with code '+
                        '{}.'.format(format))
    return out

def convert_snapshot(filename1,format1,filename2,format2,
                     overwrite=False,**kwargs):
    """
    Convert one snapshot type into another. (The old snapshot is not removed.)
        filename1: Full path to the source snapshot that should be read in.
        format1  : Integer specifying snapshot type of filename1.
        filename2: Full path to the destination snapshot that should be created.
        format2  : Integer specifying snapshot type of filename2.
        overwrite: Set to True if existing filename2 should be overwritten.
                   (default = False)
    Additional keywords are passed to the appropriate method for reading.
    """
    # Prevent overwrite
    if os.path.isfile(filename2) and not overwrite:
        raise Exception('Destination snapshot exists and overwrite is not set.')
    # Read in masses and positions from source file
    mass,pos = read_snapshot(filename1,format=format1,**kwargs)
    # Write masses and positions to the destination file
    write_snapshot(filename2,format=format2,overwrite=overwrite)
    # Return
    return

# ------------------------------------------------------------------------------
# F77 UNFORMATED BINARY SNAPSHOT FORMAT
@register_snapshot_format(0)
class UnformattedBinary(Snapshot):
    """Unformatted Fortran 77 Binary"""
    def read(self,*args,**kwargs):
        return read_unfbi77(*args,**kwargs)
    def write(self,*args,**kwargs):
        return write_unfbi77(*args,**kwargs)
    def parse_voropar(self,param):
        kwargs = {}
        if 'Unfbi77ArrayType' in param:
            kwargs['dtype'] = param['Unfbi77ArrayType']
        return kwargs

def read_unfbi77(filename,return_npart=False,dtype='float32'):
    """
    Read unformatted f77 binary snapshot.
        filename    : Full path to file that should be read.
        return_npart: If True, only the number of particles is read from the 
                      file (default = False).
        dtype       : Data type of the mass and postion arrays in the snapshot. 
                      vorovol currently only reads in float32 arrays. 
                      (default = 'float32')
    """
    import struct
    # Check data type
    dtype = np.dtype(dtype)
    if dtype!=np.dtype('float32'):
        raise Exception("vorvol assumes arrays are 32-bit floats.")
    # Open file
    fd = open(filename,'rb')
    # Read in number of particles
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != np.dtype('int32').itemsize:
        raise IOError('Error reading number of particles from file.')
    nout = struct.unpack('i',fd.read(recl))[0]
    recl = struct.unpack('i',fd.read(4))[0]
    if return_npart: return nout
    # Read in masses
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    mass = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Read in x,y,z positions
    pos = np.zeros((nout,3),dtype=np.float32)
    # X
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,0] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Y
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,1] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Z
    recl = struct.unpack('i',fd.read(4))[0]
    if recl != (nout*dtype.itemsize):
        raise IOError('Error reading x positions from file.')
    pos[:,2] = np.fromfile(fd,dtype=dtype,count=nout)
    recl = struct.unpack('i',fd.read(4))[0]
    # Close and return
    fd.close()
    return mass,pos

def write_unfbi77(filename,mass,pos,overwrite=False):
    """Write unformated f77 binary snapshot"""
    import struct
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Check data type
    dtype = np.dtype('float32')
    if mass.dtype!=dtype or pos.dtype!=dtype:
        raise Exception("vorvol assumes arrays are 32-bit floats.")
    # Open file
    fd = open(filename,'wb')
    # Write number of particles
    nout = len(mass)
    fd.write(struct.pack('i',np.dtype('int32').itemsize))
    fd.write(struct.pack('i',nout))
    fd.write(struct.pack('i',np.dtype('int32').itemsize))
    # Write masses
    fd.write(struct.pack('i',nout*dtype.itemsize))
    mass.tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # X Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,0].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Y Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,1].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Z Positions
    fd.write(struct.pack('i',nout*dtype.itemsize))
    pos[:,2].tofile(fd)
    fd.write(struct.pack('i',nout*dtype.itemsize))
    # Close and return
    fd.close()
    return

# ------------------------------------------------------------------------------
# GADGET SNAPSHOT FORMAT
# TODO
@register_snapshot_format(1)
class Gadget2(Snapshot):
    """Gadget2 Type 1 Snapshot"""
    def parse_voropar(self,param):
        kwargs = {}
        if 'GadgetParticleType' in param:
            kwargs['ptype'] = param['GadgetParticleType']
        return kwargs

# ------------------------------------------------------------------------------
# BUILDGAL SNAPSHOT FILES
# TODO?

# ------------------------------------------------------------------------------
# BUILDGAL TREEBI FILES
@register_snapshot_format(2)
class BuildgalTreebi(Snapshot):
    """Buildgal TREEBI File"""
    def read(self,*args,**kwargs):
        return read_bgtreebi(*args,**kwargs)
    def write(self,*args,**kwargs):
        return write_bgtreebi(*args,**kwargs)
    def parse_voropar(self,param):
        kwargs = {}
        if 'BgTreebiNskip' in param:
            kwargs['nskip'] = param['BgTreebiNskip']
        return kwargs

def read_bgtreebi(filename,nskip=0,return_npart=False):
    """Read Buildgal TREEBI files"""
    fd = open(filename,'r')
    # Read in header
    headline = fd.readline().strip()
    headout = headline.split()
    if headout[0] == '******':
        nout = 1000001
    else:
        nout = int(float(headout[0]))
    dim = int(float(fd.readline().strip()))
    time = float(fd.readline().strip())
    if return_npart: return nout-nskip
    # Allocate
    mass = np.zeros(nout,dtype=np.float32)
    pos = np.zeros((nout,dim),dtype=np.float32)
    # Read masses
    for i in range(nout):
        mass[i] = np.float32(fd.readline().strip())
    # Read positions
    for i in range(nout):
        posstr = fd.readline().strip().split()
        pos[i,0] = np.float32(posstr[0])
        pos[i,1] = np.float32(posstr[1])
        pos[i,2] = np.float32(posstr[2])
    # Close and return
    fd.close()
    return mass[nskip:],pos[nskip:,:]

def write_bgtreebi(filename,mass,pos,overwrite=False):
    """Write Buildgal TREEBI files"""
    # Prevent overwrite
    if os.path.isfile(filename) and not overwrite:
        print 'Specified file already exists and overwrite not set.'
        print '    '+filename
        return
    # Check sizes
    if len(mass)!=pos.shape[0]:
        raise Exception('mass has {} elements, pos has shape {}'.format(len(mass),pos.shape))
    N,dim = pos.shape
    # Open file
    fd = open(filename,'w')
    # Write header
    fd.write(' {:d}     {:d}     {:d}\n'.format(N,0,0))
    fd.write('      {:d}\n'.format(dim))
    fd.write('  {:11.6E}\n'.format(0.))
    # Write masses
    for i in range(N):
        fd.write('  {:11.6E}\n'.format(float(mass[i])))
    # Write positions
    for i in range(N):
        fd.write('  {:11.6E} {:11.6E} {:11.6E}\n'.format(float(pos[i,0]),
                                                         float(pos[i,1]),
                                                         float(pos[i,2])))
    # Close file ane return
    fd.close()
    return

