import  numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt, pow, log
from cython.parallel import prange
cdef double PI = 3.14159265359
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


def azimuthalAverage(ff):
    """
    Calculate the azimuthally averaged radial profile.
    ff - The 2D function
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

cpdef bubble(u, radi, phiP=1, phiM=-1):
    r2 = radi*radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            rsq = (i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny)
            if rsq<r2:
                u[i,j] = phiM 
            else:
                u[i,j] = phiP 
    return u


cpdef droplet(u, radi, phiP=1, phiM=-1):
    r2 = radi*radi
    Nx, Ny = np.shape(u)
    for i in range(Nx):
        for j in range(Ny):
            rsq = (i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny)
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

        
cpdef structureFactor( u, dim):
    '''
    Computes S(k) = <u(k)u(-k)> given the u(r)
    This is computed using FFT of u to obtain u(k)
    A multiplication of u(k)u(-k) is same as (abs(u(k)))^2
    if the field u is real using the definition of complex numbers
    '''
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
    

cpdef avgFunc(u, bins, dim):
    if dim==2:
        Nx, Ny = np.shape(u)
        xx, yy = np.meshgrid(np.arange(-Nx/2, Nx/2),np.arange(-Ny/2, Ny/2))
        rr = np.sqrt(xx*xx + yy*yy)
    if dim==3:
        Nx, Ny, Nz = np.shape(u)
        xx, yy, zz = np.meshgrid(np.arange(-Nx/2, Nx/2),np.arange(-Ny/2, Ny/2),np.arange(-Nz/2, Nz/2))
        rr = np.sqrt(xx*xx + yy*yy + zz*zz)
    
    rr = rr.flatten()        
    rs = np.sort(rr)
    ri = np.argsort(rr)

    u  = u.flatten();   ua = np.zeros(bins)
    u  = u[ri]

    ht, bns = np.histogram(rs, bins)
    bn = 0.5*(bns[:-1] + bns[1:])
    hm = np.cumsum(ht)

    ua[0] = np.mean( u[0:ht[0]] )
    ua[bins-1] = np.mean( u[bins-1-ht[bins-1]:] )
    for i in range(1, bins-1):
        ua[i] = np.mean( u[hm[i]+1:hm[i+1]])
    return ua, bn
