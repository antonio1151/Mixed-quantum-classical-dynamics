import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from numba import njit
from multiprocessing import Pool 
from time import time
import sys
sys.path.insert(0,"dynamics")
from Hamiltonian import hamiltonian
from EOM import EOM

#conversion paramaters
pi=np.pi
c0 = 137 #speed of light
E0 = 1/(4*pi)  #vacuum permitivity
#setting system parametrer
### Pure electronic Hamiltonian
##### Transition energies and couplings
Ener = np.diag([-0.2798,-0.6738])
Jcop = np.zeros((2,2))
###transition dipole of the electronic states
mu=1.034
mmu=np.array([[0.0,mu],[mu,0.0]])
###optical field parameters
Fn = 300 #number of modes
LC = 2.362e5 #length of the cavity
w =  np.zeros(Fn)
lam =  np.zeros(Fn)
for i in range(Fn):
    w[i] = pi*c0*(2*i+1)/LC #frequency of the modes
    lam[i] = np.sqrt(2/(E0*LC))*np.sin(pi*(2*i+1)/2)  #coupling strength
########################################
###setting up dynamics
tf =1e4 #final time in au
dt =0.01 # timestep
N =int(tf/dt) #number of steps
t = np.linspace(0,tf,N)
ns=int(7e0) #number of initial conditions for mode
n_core = 7 #number of cores
iprint=10 # number of step that are compute before being printed 
ode_q_method ='DOP853'  #integration method for quantum dynamics the options are RK45, DOP853, BDF, and LSODA
###initialization of the Hamiltonian and Dynamics functions
hh = hamiltonian(Ener,Jcop,lam,mmu,w)
dynamics = EOM()
##################################3
###initial conditions
### Quantum system: electronic part
wf0 = np.array([1+0j,0+0j],dtype=complex)
wft = np.zeros((len(wf0)),dtype=complex)
####classical system: optical field
##
t00=time()
init_x,init_p=dynamics.sampling(w,ns)
tff=time()
print("computing time of sampling:",tff-t00)
##sorting initial conditions for being use in the parallel approach
t00=time()
init_c = np.zeros((ns,Fn*2))
for i in range(ns):
    k=0
    for j in range(Fn):
        init_c[i,k] = init_x[j][i]
        init_c[i,k+1] = init_p[j][i]
        k += 2
tff=time()
print("computing time of slicing:",tff-t00)
np.savetxt("init_cond.dat",init_c)
####
''' if you already had the intial contidion for the classical system
    comment all the lines above unilt the line with "classical system" 
    and uncomment the below lines
'''
#init_c=np.loadtxt("init_cond.dat")
#init_c = init_c.reshape((ns,Fn*2))
#print(init_c.shape)
###############################################################################
#####################################
##functions to determine the quantum force
def aij(wf):
    aijj =np.real(w*lam*mu*(np.conjugate(wf[0])*wf[1]+
            np.conjugate(wf[1])*wf[0])) 
    return aijj
####################################################
##functions for multiprocessing
###dynamics
def  simulation(initial):
    wf00= wf0
    i=0
    po = np.reshape(initial,(Fn,2))
    oscilation = np.zeros((int(N/iprint)+1,2*Fn))
    quntum= np.zeros((int(N/iprint)+1,2))
    ii=0
    xname1=initial[0]
    xname2=initial[2]
    for tim in tqdm(t):
        wft[:] = wf00[:]
#        xx[:] = initial[:]
        hha =(hh.Hel()+hh.Helc(po[:,0])).flatten()
        wf = dynamics.Integ_QEOM2(wf00,dt,tim,hha,ode_q_method)
        wf00 = wf
        bij = aij(wf00)
        x = dynamics.Integ_CEOM2(po,dt,tim,bij,w,lam,mu) 
        po = np.reshape(x,(Fn,2))
        initial = np.reshape(x,(Fn*2))

        if np.mod(i,iprint) == 0:
            quntum[ii,0] = np.real(np.conjugate(wft[0])*wft[0])
            quntum[ii,1] = np.real(np.conjugate(wft[1])*wft[1])
            oscilation[ii,:] = initial[:]
            ii += 1
        i += 1
    np.save(f"dynamics_Quan_{xname1}_{xname2}",quntum)
    np.save(f"dynamics_clas_{xname1}_{xname2}",oscilation)
    del quntum,oscilation
##########################################################

pm = Pool(processes=n_core)
xxx = pm.map(simulation,init_c)
pm.close()
pm.join()

