import  numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt, pow, log
cdef double PI = 3.1415926535
from scipy.sparse import spdiags

cdef extern from "stdlib.h" nogil:
    double drand48()
    void srand48(long int seedval)

cdef extern from "time.h":
    long int time(int)
# srand48(time(0))
srand48(100)

DTYPE   = np.float
DTYP1   = np.int32
ctypedef np.float_t DTYPE_t 



cpdef structureFactor(u, dim):
    """
    * Computes S(k) = <u(k)u(-k)> given the u(r)
    * This is computed using FFT of u to obtain u(k)
    * For a real field u, the multiplication of u(k)u(-k) 
    is same as (abs(u(k)))^2


    Parameters
    ----------
    u: A field on a grid in dim dimensions
    dim: dimensionality

    Returns
    -------
        structureFactor defined (abs(u(k)))^2/N^2
    """


    if dim==1:
        uk = np.fft.fft(u)
        uk = np.fft.fftshift(uk)
        uu = np.abs(uk)

    if dim==2:
        uk = np.fft.fft2(u)
        uk = np.fft.fftshift(uk)
        uu = np.abs(uk)
    
    if dim==3:
        uk = np.fft.fftn(u)
        uk = np.fft.fftshift(uk)
        uu = np.abs(uk)

    return (uu*uu)/(np.size(u))
    



cpdef avgFuncKspace(uk, k, bins):
    """
    Obtains radial distribution of field in Fourier space
    from square of the k vector provided

    Parameters
    ----------
    uk: A field in Fourier space 
    k:  magnitude of the k field on that field 
    bins: how many bins to do

    Returns
    -------
         ua: the averaged field along radial direction 
         bn: value of the bin
    """
    
    rr = k.flatten()        
    rs = np.sort(rr)
    ri = np.argsort(rr)

    u  = uk.flatten();   ua = np.zeros(bins)
    u  = u[ri]

    ht, bns = np.histogram(rs, bins)
    bn = 0.5*(bns[:-1] + bns[1:])
    hm = np.cumsum(ht)

    ua[0] = np.mean( u[0:ht[0]] )
    ua[bins-1] = np.mean( u[bins-1-ht[bins-1]:] )
    for i in range(1, bins-1):
        ua[i] = np.mean( u[hm[i]+1:hm[i+1]])
    n1 = np.size(bn)
    return bn[1:n1-1], ua[1:n1-1]



cpdef avgFunc(u, r):
    """
    Obtains radial distribution of field 
    from Cartesian coordinates in 2D and 3D

    Parameters
    ----------
    u: A field defined on a grid of dimension dim 
    bins: how many bins in which to sum the field 
    dim: dimension

    Returns
    -------
         ua: the averaged field along radial direction 
         bn: value of the bin
    """
    
    #if dim==2:
    #    Nx, Ny = np.shape(u)
    #    xx, yy = np.meshgrid(np.arange(-Nx/2, Nx/2),np.arange(-Ny/2, Ny/2))
    #    rr = np.sqrt(xx*xx + yy*yy)
    #if dim==3:
    #    Nx, Ny, Nz = np.shape(u)
    #    xx, yy, zz = np.meshgrid(np.arange(-Nx/2, Nx/2),np.arange(-Ny/2, Ny/2),np.arange(-Nz/2, Nz/2))
    #    rr = np.sqrt(xx*xx + yy*yy + zz*zz)
    #
    #y, x = np.indices(u.shape)
    #rr = np.hypot(x, y)
    #rr = rr.flatten()        
    #rs = np.sort(rr)
    #ri = np.argsort(rr)

    #u  = u.flatten();   ua = np.zeros(bins)
    #u  = u[ri]

    #ht, bns = np.histogram(rs, bins)
    #bn = 0.5*(bns[:-1] + bns[1:])
    #hm = np.cumsum(ht)

    #ua[0] = np.mean( u[0:ht[0]] )
    #ua[bins-1] = np.mean( u[bins-1-ht[bins-1]:] )
    #for i in range(1, bins-1):
    #    ua[i] = np.mean( u[hm[i]+1:hm[i+1]])
    #return ua, bn
    #y, x = np.indices(u.shape)
    #r = np.hypot(x, y)

    #ind = np.argsort(r.flat)
    #r_sorted = r.flat[ind]
    #i_sorted = u.flat[ind]
    #r_int    = r_sorted.astype(int)

    #deltar   = r_int[1:] - r_int[:-1]  
    #rind     = np.where(deltar)[0]       
    #nr       = rind[1:] - rind[:-1]        
    #
    #csim = np.cumsum(i_sorted, dtype=float)
    #tbin = csim[rind[1:]] - csim[rind[:-1]]

    #fr = tbin / nr
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), u.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile, nr



def azimuthalAverage(ff):
    """
    Obtains radial distribution of field ff

    Parameters
    ----------
    ff: A field defined on a grid of two-dimension 

    Returns
    -------
         fr: the averaged field along radial direction 
    """
    y, x = np.indices(ff.shape)
    r = np.hypot(x, y)

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = ff.flat[ind]
    r_int    = r_sorted.astype(int)

    deltar   = r_int[1:] - r_int[:-1]  
    rind     = np.where(deltar)[0]       
    nr       = rind[1:] - rind[:-1]        
    
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    fr = tbin / nr
    return fr


def azimuthalAverageOld(ff):
    """
    Obtains radial distribution of field ff

    Parameters
    ----------
    ff: A field defined on a grid of two-dimension 

    Returns
    -------
         fr: the averaged field along radial direction 
    """
    y, x = np.indices(ff.shape)
    r = np.hypot(x, y)

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = ff.flat[ind]
    r_int    = r_sorted.astype(int)

    deltar   = r_int[1:] - r_int[:-1]  
    rind     = np.where(deltar)[0]       
    nr       = rind[1:] - rind[:-1]        
    
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    fr = tbin / nr
    return fr


def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is
             None, which then uses the center of the image (including
             fracitonal pixels).

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof


def radial_profile(data, r, bins_N=100):
    ring_brightness, radius = np.histogram(r, weights=data, bins=bins_N)
    return radius[1:], ring_brightness    




## readymade bubbles and droplets
cpdef bubble(u, radi, locx=0, locy=0, phiP=1, phiM=-1):
    r2 = radi*radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            rsq = (i-locx)*(i-locx) + (j-locy)*(j-locy)
            if rsq<r2:
                u[i,j] = phiM 
            else:
                u[i,j] = phiP 
    return u


cpdef droplet(u, radi, locx=0, locy=0, phiP=1, phiM=-1):
    r2 = radi*radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            rsq = (i-locx)*(i-locx) + (j-locy)*(j-locy)
            if rsq<r2:
                u[i,j] = phiP 
            else:
                u[i,j] = phiM 
    return u


cpdef squareDroplet(u, radi, phiP=1, phiM=-1):
    r1, r2 = radi*.6, radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            if i<(Nx/2+r1) and  i>(Nx/2-r1) and j<(Ny/2+r2) and  j>(Ny/2-r2):
                u[i,j] = phiP 
            else:
                u[i,j] = phiM 
    return u


cpdef ellipseDroplet(u, radi, phiP=1, phiM=-1):
    r1, r2 = radi*.64, radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            rsq = (i-0.5*Nx)*(i-0.5*Nx)/(r1*r1)+(j-0.5*Ny)*(j-0.5*Ny)/(r2*r2)
            if rsq<0.5:
                u[i,j] = phiP 
            else:
                u[i,j] = phiM 
    return u


cpdef twoBubbles(u, radi1, locx1, locy1, radi2, locx2, locy2, phiP=1, phiM=-1):
    rr1 = radi1*radi1
    rr2 = radi2*radi2
    Nx, Ny = np.shape(u)

    for i in range(Nx):
        for j in range(Ny):
            rsq1 = (i-locx1)*(i-locx1) + (j-locy1)*(j-locy1)
            rsq2 = (i-locx2)*(i-locx2) + (j-locy2)*(j-locy2)
            if rsq1<rr1:
                u[i,j] = phiM 
            elif rsq2<rr2:
                u[i,j] = phiM 
            else:
                u[i,j] = phiP 
    return u

##    cpdef gaussianRn(self, double fac):
##        dw = np.zeros((self.Nx, self.Ny), dtype=DTYPE)
##        cdef double [:,:] w1 = dw 
##        cdef int i, j
##
##        for i in range(self.Nx):
##            for j in range(self.Ny):
##                w1[i,j] = fac*gaussianRn()
##        return dw
##
##cdef double gaussianRn() nogil:
##    cdef int iset = 0;
##    cdef double fac, rsq, v1, v2;
##  
##    if (iset == 0): 
##        v1 = 2.0*drand48()-1.0;
##        v2 = 2.0*drand48()-1.0;
##        rsq = v1*v1 + v2*v2;
##        while (rsq >= 1.0 or rsq == 0.0):
##            v1 = 2.0*drand48()-1.0;
##            v2 = 2.0*drand48()-1.0;
##            rsq = v1*v1 + v2*v2;
##        fac = sqrt(-2.0*log(rsq)/rsq);
##        iset = 1
##        return v2*fac
##    else:
##        iset = 0
##        return v1*fac
##         
        



## finite difference using explicit loops - differs from the implmentation in dms.pyx
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class FiniteDifference:
    ''' 
    Finite Difference using loops 
    '''
    cdef:
        readonly int Np, Nx, Ny, Nz, dim, NN
        readonly double h, facy, facz  
        readonly np.ndarray iupa, jupa, idwn, jdwn

    
    def __init__(self,  grid):

        self.dim = grid.get('dim')
    
        if self.dim == 2: 
            self.Nx, self.Ny = grid.get('Nx'), grid.get('Ny')

            self.iupa  = np.empty((self.Nx), dtype=DTYP1)
            self.jupa  = np.empty((self.Ny), dtype=DTYP1)
            self.idwn  = np.empty((self.Nx), dtype=DTYP1)
            self.jdwn  = np.empty((self.Ny), dtype=DTYP1)

            for i in range(self.Nx):
                if i==self.Nx-1:
                    self.iupa[i] = 0
                else:
                    self.iupa[i] = i+1

                if i==0:
                    self.idwn[i] = self.Nx-1
                else:
                    self.idwn[i] = i-1
            
            for j in range(self.Ny):
                if j==self.Ny-1:
                    self.jupa[j] = 0
                else:
                    self.jupa[j] = j+1

                if j==0:
                    self.jdwn[j] = self.Ny-1
                else:
                    self.jdwn[j] = j-1

        elif self.dim == 3:
            self.Nx, self.Ny, self.Nz = grid.get('Nx'), grid.get('Ny'), grid.get('Nz')
                                                                         
            self.vx  = np.empty(self.Nx, self.Ny, self.Nz)
            self.vy  = np.empty(self.Nx, self.Ny, self.Nz)
            self.vz  = np.empty(self.Nx, self.Ny, self.Nz)

    cpdef diffx(self, double [:,:] u):        
        '''
        computes d/dx u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j,
        cdef double [:,:] du1 = du 
        cdef int [:] iup = self.iupa
        cdef int [:] jup = self.jupa
        cdef int [:] idw = self.idwn
        cdef int [:] jdw = self.jdwn 

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = .1*(u[iup[i],jup[j]]-u[iup[i],jdw[j]]) + .3*(u[i,jup[j]]-u[i,jdw[j]]) +.1*(u[idw[i],jup[j]]-u[idw[i],jdw[j]])
        return du

    cpdef diffy(self, double [:,:] u):        
        '''
        computes d/dy u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j,
        cdef double [:,:] du1 = du 
        cdef int [:] iup = self.iupa
        cdef int [:] jup = self.jupa
        cdef int [:] idw = self.idwn
        cdef int [:] jdw = self.jdwn 

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = .1*(u[iup[i],jup[j]]-u[idw[i],jup[j]]) + .3*(u[iup[i],j]-u[idw[i],j]) +.1*(u[iup[i],jdw[j]]-u[idw[i],jdw[j]])
        return du

    cpdef laplacian(self, double [:,:] u):        
        '''
        computes D^2 u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j
        cdef double [:,:] du1 = du 
        cdef int [:] iup = self.iupa
        cdef int [:] jup = self.jupa
        cdef int [:] idw = self.idwn
        cdef int [:] jdw = self.jdwn

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = -0.5*u[idw[i],jup[j]]+2.0*u[i,jup[j]]-0.5*u[iup[i],jup[j]] + 2.0*u[idw[i],j]-6.0*u[i,j]+2.0*u[iup[i],j] - 0.5*u[idw[i],jdw[j]]+2.0*u[i,jdw[j]]-0.5*u[iup[i],jdw[j]]
        return du
    
    cpdef diffx1(self, double [:,:] u):        
        '''
        computes d/dx u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j, jup, jup2, jup3, jup4, jup5, jdw, jdw2, jdw3, jdw4, jdw5
        cdef double [:,:] du1 = du 
        cdef int [:] iupa = self.iupa
        cdef int [:] jupa = self.jupa
        cdef int [:] idwn = self.idwn
        cdef int [:] jdwn = self.jdwn 

        cdef double fp = 4.0/105, tt = 1.0/280.0

        for i in range(self.Nx):
            for j in range(self.Ny):
                jup  = jupa[j]
                jup2 = jupa[jupa[j]]
                jup3 = jupa[jupa[jupa[j]]]
                jup4 = jupa[jupa[jupa[jupa[j]]]]
                jup5 = jupa[jupa[jupa[jupa[iupa[j]]]]]
                jdw  = jdwn[j];
                jdw2 = jdwn[jdwn[j]];
                jdw3 = jdwn[jdwn[jdwn[j]]];
                jdw4 = jdwn[jdwn[jdwn[jdwn[j]]]];
                jdw5 = jdwn[jdwn[jdwn[jdwn[jdwn[j]]]]];

                du1[i,j] = -tt*u[i,jup4]+fp*u[i,jup3]-.2*u[i,jup2]+.8*u[i,jup]+tt*u[i,jdw4]-fp*u[i,jdw3]+.2*u[i,jdw2]-.8*u[i,jdw]
        return du 

    cpdef diffx2(self, double [:,:] u):        
        '''
        computes d/dx u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j, jup, jup2, jup3, jup4, jup5, jdw, jdw2, jdw3, jdw4, jdw5
        cdef double [:,:] du1 = du 

        cdef double t1 = 1.0/3, t2 = 1.0/12

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = t1*(u[i,j+1]-u[i,j-1]) + t2*(u[i+1,j+1]+u[i+1,j-1]) - t2*(u[i-1,j+1]+u[i-1,j-1])
        return du 
    
    cpdef diffy2(self, double [:,:] u):        
        '''
        computes d/dx u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j, jup, jup2, jup3, jup4, jup5, jdw, jdw2, jdw3, jdw4, jdw5
        cdef double [:,:] du1 = du 

        cdef double t1 = 1.0/3, t2 = 1.0/12

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = t1*(u[i+1,j]-u[i-1,j]) + t2*(u[i+1,j-1]+u[i+1,j+1]) - t2*(u[i-1,j-1]+u[i-1,j+1])
        return du 
    
    cpdef lap2(self, double [:,:] u):        
        '''
        computes d/dx u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j, jup, jup2, jup3, jup4, jup5, jdw, jdw2, jdw3, jdw4, jdw5
        cdef double [:,:] du1 = du 

        cdef double t1 = 1.0/6, t2 = 4.0/6, t3=20/6.0

        for i in range(self.Nx):
            for j in range(self.Ny):
                du1[i,j] = t1*(u[i-1,j-1]+u[i-1,j+1]+u[i+1,j+1]+u[i+1,j-1]) + t2*( u[i,j+1]+u[i,j-1]+u[i-1,j]+u[i+1,j] ) - t3*u[i,j]
        return du 

    cpdef diffy1(self, double [:,:] u):        
        '''
        computes d/dy u
        '''
        du = np.zeros(np.shape(u), dtype=DTYPE)
        cdef int i, j, iup, iup2, iup3, iup4, iup5, idw, idw2, idw3, idw4, idw5
        cdef double [:,:] du1 = du 
        cdef int [:] iupa = self.iupa
        cdef int [:] jupa = self.jupa
        cdef int [:] idwn = self.idwn
        cdef int [:] jdwn = self.jdwn 

        cdef double fp = 4.0/105, tt = 1.0/280.0

        for i in range(self.Nx):
            iup  = iupa[i]
            iup2 = iupa[iupa[i]]
            iup3 = iupa[iupa[iupa[i]]]
            iup4 = iupa[iupa[iupa[iupa[i]]]]
            iup5 = iupa[iupa[iupa[iupa[iupa[i]]]]]
            idw  = idwn[i];
            idw2 = idwn[idwn[i]];
            idw3 = idwn[idwn[idwn[i]]];
            idw4 = idwn[idwn[idwn[idwn[i]]]];
            idw5 = idwn[idwn[idwn[idwn[idwn[i]]]]];
            for j in range(self.Ny):
                du1[i,j] =  -tt*u[iup4,j]+fp*u[iup3,j]-.2*u[iup2,j]+.8*u[iup,j]+tt*u[idw4,j]-fp*u[idw3,j]+.2*u[idw2,j]-.8*u[idw,j]
        return du
