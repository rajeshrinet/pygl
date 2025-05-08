from scipy.io import savemat
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from tqdm import tqdm
from matplotlib.animation import FuncAnimation

import scipy.fft
FFT = scipy.fft.fft2
IFFT = scipy.fft.ifft2

fft2  = np.fft.fft2
ifft2 = np.fft.ifft2
randn = np.random.randn
    
    
    
class FPS_ModelB():
    '''Class to simulate field theories using a PSEUDO-SPECTRAL TECHNIQUE WITH FILTERING'''
    def __init__(self, param):
        self.Nt = param['Nt'];    self.dt=param['dt'];    self.Nf = param['Nf']
        self.h  = param['h'];     self.a = param['a'];    self.b  = param['b']
        self.kp = param['kp'];    self.k_tilde= param['k_tilde'];  self.Ng= param['Ng']
        self.Df = 1j*param['Df']; Ng = self.Ng
        self.k_c= param['k_c']
        
        self.XX  = np.zeros((int(self.Nf+1), Ng*Ng)); 
        self.LL  = np.zeros((int(self.Nf+1), Ng*Ng)); 
        
        qf=np.fft.fftfreq(Ng)*(2*np.pi/self.h)
        self.qx, self.qy = np.meshgrid(qf, qf) 
        self.Dx, self.Dy = 1j*self.qx, 1j*self.qy

        self.q2 = self.qx*self.qx + self.qy*self.qy;    
        alpha_q = -self.q2*(self.a + self.kp*self.q2)*self.dt;   
        self.eA = np.exp(alpha_q);   self.q2b=self.q2*self.b 
        
        iq2=1/(self.q2+1e-16);  iq2[0,0]=0; self.iq2=iq2   
        self.iQ2x, self.iQ2y = self.qx*iq2, self.qy*iq2
        self.qmod = np.sqrt( (self.q2) )

        freqs = (2*np.pi)*np.fft.fftfreq(Ng)
        self.kx_s = np.fft.fftshift(freqs)
        
    def integrate(self, u):
        '''  simulates the equation and plots it at different instants '''
        ii=0;  t=0;  dt=self.dt;    Df=self.Df; simC=(100/self.Nt)
        self.u=u
        for i in tqdm(range(self.Nt)):          
            self.rhs()
            if i%(int(self.Nt/self.Nf))==0:  
                self.XX[ii,:] = (np.real(ifft2(self.u))).flatten()
                ii += 1   
            
        
    def rhs(self):
        '''
        integrator to simulate the dynamics of order paramter 
        using filtere-pseudo-spectral (FPS)
        '''       
        ## evolve the order paramter using FPS
        uc  = ifft2(self.u);    
        N_u = -self.q2b*fft2(uc*uc*uc)       

        self.u = (self.u + N_u*self.dt)*(self.eA) 
        return        

    
    def getLength(self):
        u1 = np.abs( (self.u) );  
        a2, b2 = avgFuncKspace(u1*u1/(Ng*Ng), self.qmod, int(self.Ng/4))
        
        LL = 2*np.pi*np.sum(b2)/np.sum(a2*b2)
        print (np.max(a2), np.max(b2), LL)
        return LL


    def structureFactor(self, u, dim):
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
    
        return (uu*uu) 
    
    
    
    def avgFunc(self, u, bins, dim):
        if dim==2:
            Nx, Ny = np.shape(u)
            # xx, yy = np.meshgrid(np.arange(-Nx/2, Nx/2),np.arange(-Ny/2, Ny/2))
            # rr = np.sqrt(xx*xx + yy*yy)
            kx = self.kx_s
            kx, ky = np.meshgrid(kx,kx);
            k2  = kx*kx + ky*ky; k4=k2*k2 ; 
            ksq= np.fft.fftshift((k2))
            rr = np.sqrt(ksq)
        
        if dim==3:
            Nx, Ny, Nz = np.shape(u)
            kx = self.kx_s
            xx, yy, zz = np.meshgrid(kx, kx, kx)
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

