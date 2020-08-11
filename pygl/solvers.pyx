import  numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
cdef double PI = 3.14159265359


DTYPE   = np.float
ctypedef np.float_t DTYPE_t
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef class Stokes:
    """
    Numerical solution of Stokes equation on 2D or 3D grid. 

    ...

    Parameters
    ----------
    eta: float 
        Viscosity of the fluid (eta)
    grid: dict 
        Grid and its properties as a dict

    Example
    ----------
    >>> import pygl 
    >>> eta    = .1
    >>> grid   = {"dim":2, "Nx":32, "Ny":32}
    >>> stokes = pygl.solvers.Stokes(eta, grid)

   """

    
    cdef:
        readonly int Np, Nx, Ny, Nz, dim, NN
        readonly double facx, facy, facz , eta
        readonly np.ndarray vx, vy, vz, fkx, fky, fkz, vkx, vky, vkz
    
    def __init__(self, eta, grid):

        self.dim = grid.get('dim')
        self.eta = eta
    
        if self.dim == 2: 
            self.Nx, self.Ny = grid.get('Nx'), grid.get('Ny')
            self.facx = 2*PI/self.Nx
            self.facy = 2*PI/self.Ny

            self.vx  = np.empty((self.Nx, self.Ny))
            self.vy  = np.empty((self.Nx, self.Ny))
            self.fkx = np.empty((self.Nx, self.Ny), dtype=np.complex128)    
            self.fky = np.empty((self.Nx, self.Ny), dtype=np.complex128)
            self.vkx = np.empty((self.Nx, self.Ny), dtype=np.complex128)    
            self.vky = np.empty((self.Nx, self.Ny), dtype=np.complex128)
        
        elif self.dim == 3:
            self.Nx, self.Ny, self.Nz = grid.get('Nx'), grid.get('Ny'), grid.get('Nz')
            self.facx = 2*PI/self.Nx
            self.facy = 2*PI/self.Ny
            self.facz = 2*PI/self.Nz
                                                                         
            self.vx  = np.empty((self.Nx, self.Ny, self.Nz))
            self.vy  = np.empty((self.Nx, self.Ny, self.Nz))
            self.vz  = np.empty((self.Nx, self.Ny, self.Nz))
            self.fkx = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128) 
            self.fky = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128)
            self.fkz = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128)
            self.vkx = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128) 
            self.vky = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128)
            self.vkz = np.empty((self.Nx, self.Ny, self.Nz), dtype=np.complex128)
       
        else:
            raise Exception('Current support is only for 2D or 3D')
       
       
         
    cpdef solve(self, fk):
        """
        Compute flow given force per unit area 

        self.vx, self.vy and self.vz contains the flow computed

        ...

        Parameters
        ----------
        fk: np.array 
            Fourier transform of the force per unit area on the grid

        Example
        ----------
        >>> import pygl 
        >>> eta    = .1
        >>> grid   = {"dim":2, "Nx":32, "Ny":32}
        >>> stokes = pygl.solvers.Stokes(eta, grid)
        >>> fkx    = np.random.random((32, 32))
        >>> fky    = np.random.random((32, 32))
        >>> stokes.solve( (fkx, fky) )
       """

        cdef int dim = self.dim
        if dim == 2:
            self._solve2d(fk)

        elif dim == 3:
            self._solve3d(fk)
        return
        
        
    cdef _solve2d(self, fk):        
         
        self.fkx = fk[0]
        self.fky = fk[1]        
        
        cdef:
            complex [:,:] fkx = self.fkx
            complex [:,:] fky = self.fky
            complex [:,:] vkx = self.vkx
            complex [:,:] vky = self.vky
            
            int jx, jy, Nx, Ny
            double facx, facy, iksq, kx, ky, ieta
            double complex fdotk
        Nx = self.Nx
        Ny = self.Ny 
    
        facx = self.facx
        facy = self.facy
        ieta = 1.0/self.eta
        
        for jy in range(Ny): 
            ky = jy*facy if jy<=Ny/2 else (-Ny+jy)*facy
            for jx in range(Nx): 
                kx = jx*facx if jx<=Nx/2 else (-Nx+jx)*facx
                if kx == 0 and ky == 0:
                    iksq = 0.0
                else:
                    iksq = 1.0/(kx*kx + ky*ky)
                fdotk    = kx*fkx[jy, jx] + ky*fky[jy, jx]
                
                vkx[jy, jx] = ( fkx[jy, jx] - fdotk*kx*iksq )*iksq*ieta
                vky[jy, jx] = ( fky[jy, jx] - fdotk*ky*iksq )*iksq*ieta

        vkx[0, 0] = 0.0     #the integral of the fluid flow vanishes
        vky[0, 0] = 0.0
        
        self.vx = np.real(np.fft.ifft2(self.vkx))
        self.vy = np.real(np.fft.ifft2(self.vky))
        return


    cdef _solve3d(self, fk):
        self.fkx = fk[0]
        self.fky = fk[1]        
        self.fkz = fk[2]        
        
        cdef:
            complex [:,:,:] fkx = self.fkx
            complex [:,:,:] fky = self.fky
            complex [:,:,:] fkz = self.fkz
            complex [:,:,:] vkx = self.vkx
            complex [:,:,:] vky = self.vky
            complex [:,:,:] vkz = self.vkz
            
            int jx, jy, jz, Nx, Ny, Nz
            double facx, facy, facz, iksq, kx, ky, kz, ieta
            double complex fdotk
        Nx = self.Nx
        Ny = self.Ny 
        Nz = self.Nz 
    
        facx = self.facx
        facy = self.facy
        facz = self.facz
        ieta = 1.0/self.eta
        
        for jz in range(Nz): 
            for jy in range(0, Ny): 
                for jx in range(0, Nx): 
                    kx = jx*facx if jx <= Nx / 2 else (-Nx+jx)*facx
                    ky = jy*facy if jy <= Ny / 2 else (-Ny+jy)*facy
                    kz = jz*facz if jz <= Nz / 2 else (-Nz+jz)*facz

                    if kx == 0 and ky == 0:
                        iksq = 0.0
                    else:
                        iksq = 1.0/(kx*kx + ky*ky)
                    fdotk = kx*fkx[jz, jy, jx] + ky*fky[jz, jy, jx] + kz*fkz[jz, jy, jx]
                    
                    vkx[jy, jx, jz] = ( fkx[jy, jx, jz] - fdotk*kx*iksq )*iksq*ieta
                    vky[jy, jx, jz] = ( fky[jy, jx, jz] - fdotk*ky*iksq )*iksq*ieta 
                    vkz[jy, jx, jz] = ( fkz[jy, jx, jz] - fdotk*kz*iksq )*iksq*ieta 

        vkx[0, 0, 0] = 0.0  #the integral of the fluid flow vanishes
        vky[0, 0, 0] = 0.0
        vkz[0, 0, 0] = 0.0
        
        self.vx = np.real(np.fft.ifftn(self.vkx))
        self.vy = np.real(np.fft.ifftn(self.vky))
        self.vz = np.real(np.fft.ifftn(self.vkz))

        return

