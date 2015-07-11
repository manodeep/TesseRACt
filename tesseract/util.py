import os
import numpy as np
import sys

# ------------------------------------------------------------------------------
# FITTING UTILITIES
def myleastsq(errfunc0,x0,args=None,bounds=None,**exkw):
    """Hack to allow leastsq to take bounds if minimize function is missing"""
    from scipy import optimize
    if hasattr(optimize,'minimize'):
        def errfunc(x,*iargs):
            return sum(errfunc0(x,*iargs)**2)
        if args is not None: exkw['args'] = args
        res = optimize.minimize(errfunc,x0[:],bounds=bounds,**exkw)
        return res.x,res.success
    else:
        lres = sys.float_info.max
        def errfunc(x,*iargs):
            if bounds!=None:
                for idx in range(len(x)):
                    if bounds[idx][0]!=None and x[idx]<bounds[idx][0]: return lres
                    if bounds[idx][1]!=None and x[idx]>bounds[idx][1]: return lres
            return errfunc0(x,*iargs)
        if args is not None: exkw['args'] = args
        return optimize.leastsq(errfunc,x0,**exkw)

# ------------------------------------------------------------------------------
# UTILITIES FOR CALCULATING SPHERICAL THINGS
def pos2rad(pos):
    """Calculates radius from position"""
    return np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
def sphvol(r):
    """Returns spherical volume for provided radii"""
    return (4./3.)*np.pi*(r**3.)
def sphrad(vol):
    """Returns spherical radii for provided volumes"""
    return (3.*vol/(4.*np.pi))**(1./3.)
def calc_vcirc(r,menc,G=1.):
    """Returns the circular velocity"""
    if G is None: G = 1.
    return np.sqrt(G*menc/r)
def calc_menc(m,sortby=None):
    """
    Calculates the total mass enclosed by particles from their radii and mass.
        m     : masses of particles
        sortby: quantity that masses should be sorted by before summing. The
                enclosed mass will be returned in the original particle order.
    """
    # Sort and sum if sortby array provided
    if sortby is not None:
        idxsort = np.argsort(sortby)
        idxsort_rev = np.argsort(idxsort)
        menc = np.cumsum(m[idxsort])[idxsort_rev]
    # Otherwise just sum
    else:
        menc = np.cumsum(m)
    # Return
    return menc

def calc_rhoenc(mass,r,rmax):
    """Returns the density enclosed within rmax"""
    idx = (r<rmax)
    return mass[idx].sum()/sphvol(rmax)

def vol2rad(vol,sortby=False,outsort=False,weightby=False,weight=0.5):
    """
    Calculate particle radii from volume information. Volumes are sorted
    then cummulatively summed to get enclosed volume at each particle. Radii are
    then calculated by assumming that each enclosed volume is a sphere. Returned
    radii are in the same particle order as the input volumes.
        vol     : (N,) array of volumes corresponding to each particle
        sortby  : (N,) array to sort volumes by before calculating radii. If not
                  provided, vol is used. If set to True, volumes are assumed to 
                  be already sorted.
        outsort : if True, the array for sorting the particles is also output
                  and the returned radii are already sorted.
        weightby: (N,) array of values to weight volumes by after calculating
                  radii. Sorting occurs before weighting. This is useful if
                  you would like to weight by the real radius or another 
                  measure of radius. The returned radii are then the weighted
                  arithmetic mean of the volume based radii (rvol) and the
                  weightby array (i.e. weight*weightby + rvol*(1-weightby))
                  The resulting array is again sorted if outsort is True.
                  (default = False)
        weight  : float value used to weight the weightby array by. Only used
                  if weightby is not False. (default = 0.5)
    """
    # Sort
    if isinstance(sortby,np.ndarray):
        if len(sortby)!=len(vol):
            raise Exception('The sortby array (len={}) must have the '.format(len(sortby))+
                            'same length as the volume array (len={}).'.format(len(vol)))
        idx = np.argsort(sortby)
    elif isinstance(sortby,bool):
        if sortby:
            idx = slice(0,len(vol))
        else:
            idx = np.argsort(vol)
    else:
        raise Exception('Invalid sortby value: {}'.format(sortby))
    # Sort volumes and sum
    rad = sphrad(np.cumsum(vol[idx]))
    # Weight radii
    if isinstance(weightby,np.ndarray):
        if len(weightby)!=len(vol):
            raise Exception('The weightby array (len={}) must have the '.format(len(weightby))+
                            'same length as the volume array (len={}).'.format(len(vol)))
        rad = weight*weightby + (1.-weight)*rad[np.argsort(idx)]
        if outsort:
            idx = np.argsort(rad)
            rad = rad[idx]
            return rad,idx
        else:
            return rad
    # Return
    if outsort:
        return rad,idx
    else:
        if sortby == True:
            return rad
        else:
            idxrev = np.argsort(idx)
            return rad[idxrev]

# ------------------------------------------------------------------------------
# TRANSFORMATION THINGS
def squeeze_constvol(pos,yfact=-1,zfact=-1):
    """
    Squeezes 3D positions along y and z axes while preserving volume.
        pos  : (N,3) array of positions
        yfact: factor that y direction is squeezed by (0 < yfact < 1)
        zfact: factor that z direction is squeezed by (0 < zfact < 1)
    """
    squfact = 1.0
    if yfact > 0:
        squfact*=yfact
        pos[:,1]*=yfact
    if zfact > 0:
        squfact*=zfact
        pos[:,2]*=zfact
    pos*=(squfact**(-1./3.))
    return pos

def triax2squeeze(t,e):
    """
    Return y and z squeeze factors for given triaxiality and ellipticity 
        t: triaxiality defined as (c**2 - b**2)/(c**2 - a**2) where a > b > c
        e: overall ellipticity defined as c/a
    """
    a = 1.
    c = e
    b = np.sqrt(c**2 - t*(c**2-a**2))
    return b

# ------------------------------------------------------------------------------
# STRING THINGS
def val2str(val):
    """Writes values to a string"""
    # Return the input if it's a string
    if   isinstance(val,str  ): valstr=val
    # Handle types where spaces are added
    elif isinstance(val,tuple): valstr=repr(val).replace(', ',',')
    elif isinstance(val,list ): valstr=repr(val).replace(', ',',')
    elif isinstance(val,dict ): valstr=repr(val).replace(', ',',').replace(': ',':')
    # Otherwise use repr()
    else: valstr=repr(val)
    # Return output
    return valstr

def str2val(valstr,dtype=None):
    """Recovers variables from a string"""
    valstr_strip=valstr.strip()
    if isinstance(dtype,type): dtype = str(dtype)
    # Get data type if not provided
    if not isinstance(dtype,str):
        if len(valstr_strip)==0: dtype='str'
        elif valstr_strip.startswith('"') and valstr_strip.endswith('"'):
            valstr_strip=valstr_strip.strip('"')
            dtype='str'
        elif valstr_strip.lower() == 'nan': dtype='float'
        elif valstr_strip.endswith('L') and valstr_strip.rstrip('L').isdigit(): 
            dtype='long'
        elif valstr_strip.isdigit() or (valstr_strip[0] in ['-','+'] and valstr_strip[1:].isdigit()):
            if int(valstr_strip) == long(valstr_strip): dtype='int'
            else                                      : dtype='long'
        elif valstr_strip.startswith("'") and valstr_strip.endswith("'"):
            valstr_strip=valstr_strip.strip("'")
            dtype='str'
        elif valstr_strip.lower() == 'none'          : dtype='none'
        elif valstr_strip.strip('.').lower() in ['true','false']: dtype='bool'
        elif valstr_strip.startswith('(') and valstr_strip.endswith(')'):
            if ',' in valstr: dtype='tuple'
            elif 'j' in valstr: dtype='complex'
            else: dtype='str'
        elif valstr_strip.startswith('[') and valstr_strip.endswith(']'): 
            dtype='list'
        elif valstr_strip.startswith('{') and valstr_strip.endswith('}'): 
            dtype='dict'
    # If you know the data type, convert it
    if isinstance(dtype,str):
        dtype=dtype.lower()
        if   dtype=='s' or 'str' in dtype:
            if   valstr_strip.startswith("'") and valstr_strip.endswith("'"): 
                val=valstr_strip.strip("'")
            elif valstr_strip.startswith('"') and valstr_strip.endswith('"'): 
                val=valstr_strip.strip('"')
            else: val=valstr_strip
        elif dtype=='f' or 'float'   in dtype: val=float(valstr_strip)
        elif dtype=='d' or 'double'  in dtype: val=float(valstr_strip)
        elif dtype=='i' or 'int'     in dtype: 
            val=int(float(valstr_strip.replace('L','')))
        elif dtype=='l' or 'long'    in dtype: 
            val=long(float(valstr_strip.replace('L','')))
        elif dtype=='c' or 'complex' in dtype: 
            val=complex(valstr_strip.replace(' ',''))
        elif dtype=='n' or 'none'    in dtype: val=None
        elif dtype=='b' or 'bool'    in dtype:
            if   'true'  in valstr_strip.strip('.').lower(): val=True
            elif 'false' in valstr_strip.strip('.').lower(): val=False
            else: 
                raise Exception('Unrecognised string for bool type: {}'.format(valstr))
        elif 'tuple' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val=[]
            for ivalstr in vallist:
                val.append(str2val(ivalstr))
            val=tuple(val)
        elif 'list' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val=[]
            for ivalstr in vallist:
                val.append(str2val(ivalstr))
        elif 'dict' in dtype:
            vallist=valstr_strip[1:-1].split(',')
            val={} ; val0={} ; key0=[]
            for ivalstr in vallist:
                if ':' in ivalstr:
                    ikey,ival = ivalstr.split(':')
                    val0[ikey]=ival
                    key0.append(ikey)
                else:
                    val0[key0[-1]]+=','+ivalstr
    # If you don't know the data type try a few things
    else:
        try:
            # Try as complex
            if 'j' in valstr_strip:
                val=complex(valstr_strip.replace(' ',''))
            # Try as float
            else:
                val=float(valstr_strip)
        # Otherwise return as string
        except ValueError:
            val=valstr_strip
    # Return output
    return val

def num2str(numval,formstr=None,decrep='p',negrep='n',nsig=3,trimzero=True):
    """
    Convert number to string, replacing punctuation.
        numval  : value to be converted
        formstr : format string used for the conversion
        decrep  : string that '.'s are replaced with (default = 'p')
        negrep  : string that '-'s are replaced with (default = 'n')
        nsig    : number of significan digits (default = 3)
        trimzero: if true, remove trailing zeros in decimal (default = True)
    """
    import copy
    from string import maketrans
    # Set default
    if formstr is None: formstr='.{}g'.format(nsig)
    # Round numval if within precision
    if trimzero and not isinstance(numval,(int,long)):
        precval=10.**(-(nsig-1))
        if numval == 0:
            numval=int(numval)
        else:
            if int(numval) == 0:
                if abs(numval) < precval:
                    numval=int(numval)
            else:
                if abs(numval % int(numval)) < precval: 
                    numval=int(numval)
    # Set default format for integers/longs
    if isinstance(numval,(int,long)): formstr='d'
    # Create string
    decstr=(r'{:'+formstr+'}').format(numval)
    # Trim zeros
    if trimzero and '.' in decstr:
        if 'e' in decstr or 'E' in decstr:
            if 'e' in decstr: decbas,decexp=decstr.split('e')
            if 'E' in decstr: decbas,decexp=decstr.split('E')
            if 'e' and 'E' in decstr:
                raise Exception('[mmlstring.dec2str] Both e and E are in the string: {}'.format(decstr))
        else:
            decbas=copy.deepcopy(decstr)
            decexp=''
        decbas.rstrip('0')
        decbas.rstrip('.')
        decstr=decbas+decexp
    # Create tables to replace .- characters and remove + .
    rtab='+ '
    itab='.-'
    ftab=decrep+negrep
    trantab=maketrans(itab,ftab)
    # Return translated string
    return decstr.translate(trantab,rtab)


