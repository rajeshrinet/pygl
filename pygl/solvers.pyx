import  numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt, pow, log, sin, cos, exp
from cython.parallel import prange

cdef double PI = 3.14159265359
fft2  = np.fft.fft2
ifft2 = np.fft.ifft2
randn = np.random.randn

cdef extern from "stdlib.h" nogil:
    double drand48()
    void srand48(long int seedval)

cdef extern from "time.h":
    long int time(int)
srand48(time(0))

cdef double gaussianRn() nogil:
    cdef int iset = 0;
    cdef double fac, rsq, v1, v2;
  
    if (iset == 0): 
        v1 = 2.0*drand48()-1.0;
        v2 = 2.0*drand48()-1.0;
        rsq = v1*v1 + v2*v2;
        while (rsq >= 1.0 or rsq == 0.0):
            v1 = 2.0*drand48()-1.0;
            v2 = 2.0*drand48()-1.0;
            rsq = v1*v1 + v2*v2;
        fac = sqrt(-2.0*log(rsq)/rsq);
        iset = 1
        return v2*fac
    else:
        iset = 0
        return v1*fac


def getRN():
    return (gaussianRn())

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



@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef class ModelB():
    """
    Simualting model B
    """
    cdef readonly int Nt, Ny, Nz, Ng, Nf, Tf
    cdef readonly np.ndarray kx, ky, ksq, XX, du, alp
    cdef readonly double dt, h, a, b, kp, Df, t

    def __init__(self, param):

        self.Nt = param['Nt']
        self.dt = param['dt']
        self.Nf = param['Nf']
        self.h  = param['h']
        self.a  = param['a']
        self.b  = param['b']
        self.t  = 0.0
        self.kp = param['kp']

        self.Ng = param['Ng']
        self.Tf =int(self.Nt/self.Nf)
        self.Df = param['Df']
        Ng = self.Ng
        
        self.XX  = np.zeros((int(self.Nf), Ng*Ng)) 
        self.du  = np.zeros((Ng,Ng), dtype=np.complex128) 
        
        kx  = np.fft.fftfreq(Ng)*(2*np.pi/self.h)
        ky  = np.fft.fftfreq(Ng)*(2*np.pi/self.h)
        self.kx, self.ky = np.meshgrid(kx, ky) 
        self.ksq = self.kx*self.kx + self.ky*self.ky 
        self.alp = np.exp(-self.dt*self.ksq*(self.a + self.kp*self.ksq))
    
    
    cpdef integrate(self, u):
        '''  simulates the equation and plots it at different instants '''
        cdef int ii=0, i ,Nt=self.Nt, Tf=self.Tf,
        cdef double  simC=(100.0/self.Nt)

        for i in range(Nt):          
            self.rhs(u)
            u = u + self.du

            if i%(Tf)==0:  
                self.XX[ii,:]=(np.real(ifft2(u))).flatten()
                ii += 1   
                if ii%50==0:
                    print (int(i*simC), '% done')
 
     
    cpdef integrateFPS(self, u):
        '''  simulates the equation and plots it at different instants '''
        cdef int ii=0, i ,Nt=self.Nt, Tf=self.Tf, Ng=self.Ng
        cdef double  simC=(100.0/self.Nt), dt=self.dt, t
        
        self.du = u
        for i in range(Nt):          
            self.rhsFPS()

            if i%(Tf)==0:  
                self.XX[ii,:]=(np.real(ifft2(self.du))).flatten()
                ii += 1   
                if ii%50==0:
                    print (int(i*simC), '% done')
                
    cpdef rhs(self, uk):
        '''
        returns the right hand side of \dot{phi} in model B
        \dot{phi} = Δ(a*u + b*u*u*u + kΔu + λ(∇u)^2) 
        '''
        cdef int i, j, Ng=self.Ng 
        cdef double a=self.a, b=self.b, kp=self.kp, k2, dt=self.dt
        cdef double [:, :] kx = self.kx
        cdef double [:, :] ky = self.ky
        cdef double [:, :] ksq = self.ksq
        cdef complex [:, :] du  = self.du
        cdef complex [:, :] u   = uk 
        cdef complex Df = 1j*self.Df
        
        cdef double [:, :] rnx = randn(Ng, Ng)
        cdef double [:, :] rny = randn(Ng, Ng)
        
        uc =ifft2(uk)
        cdef complex [:,:] u3 = (fft2(uc*uc*uc))

        for i in prange(Ng, nogil=True):
            for j in range(Ng):
                k2 = ksq[i,j]
                du[i,j] = -k2*dt*( (a+kp*k2)*u[i,j] + b*u3[i,j] ) + Df*(kx[i,j]*rnx[i,j]+ky[i,j]*rny[i,j])
        return 
    
    
    cpdef rhsFPS(self):
        '''
        returns the right hand side of \dot{phi} in model B
        \dot{phi} = Δ(a*u + b*u*u*u + kΔu + λ(∇u)^2) 
        '''
        cdef int i, j, Ng=self.Ng 
        cdef double a=self.a, kp=self.kp
        cdef double kx, ky, k2, dtb=self.b*self.dt
        cdef complex Df = 1j*self.Df

        cdef double [:, :] kxv = self.kx
        cdef double [:, :] kyv = self.ky
        cdef double [:, :] rnx = randn(Ng, Ng)
        cdef double [:, :] rny = randn(Ng, Ng) 
        cdef double [:, :] alp = self.alp
        cdef complex [:, :] du = self.du 

        
        uc =ifft2(self.du) 
        cdef complex [:,:] u3 = (fft2(uc*uc*uc))

        for i in prange(Ng, nogil=True):
            for j in range(Ng):
                kx = kxv[i,j]
                ky = kyv[i,j]
                k2 = kx*kx + ky*ky
                du[i,j] = (du[i,j] - k2*dtb*u3[i,j] + Df*(kx*rnx[i,j]+ky*rny[i,j]))*(alp[i,j])
        return 
