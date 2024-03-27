import numpy as np
from scipy.integrate import odeint,ode,solve_ivp

class EOM(object):
    def __init__(self):
        super().__init__()
        self.hbar = 1.0#0.6582 #plank constant
        self.kb = 1.0 #8.617333E-5 #Boltzmann constant
        self.time_nor_nat=0.6582 #from natural units to fs
#Quantum equation of motion 
    def QEOM(self,dt,wf,*args):
        Hamil = np.reshape(args,(len(wf),len(wf)))
        return np.dot(-1j*Hamil/self.hbar,wf)
    def Integ_QEOM(self,wf,dt,t,*args):
        k1 = dt*self.QEOM(t,wf,args)
        k2 = dt*self.QEOM(t+dt/2,wf+k1/2,args)
        k3 = dt*self.QEOM(t+dt/2,wf+k2/2,args)
        k4 = dt*self.QEOM(t+dt,wf+k3,args)
        y = wf+(1/6)*(k1+2*(k2+k3)+k4)
        return y        
    def Integ_QEOM2(self,wf,dt,t,Hamil,ODE_method):
        yt = solve_ivp(self.QEOM,[t,t+dt],wf,method=ODE_method,args=Hamil)
        return yt.y[:,-1]

# the Classical equation of motion
##aij is the quantum component of the feedback force
    def CEOM(self,t,y,aij,w,lam,mu):
        v = y[:,1]
        Ac = -w**2*y[:,0]-w*lam*mu*aij
        return np.column_stack((v,Ac))
    def CEOM2(self,t,y,aij,w,lam,mu):
        y = np.reshape(y,(len(w),2))
        v = y[:,1]
        Ac = -w**2*y[:,0]-w*lam*mu*aij
        return np.column_stack((v,Ac)).flatten()
# the integrator is the RK4
    def Integ_CEOM(self,y0,dt,t,aij,w,lam,mu):
        k1 = dt*self.CEOM(t,y0,aij,w,lam,mu)
        k2 = dt*self.CEOM(t+dt/2,y0+k1/2,aij,w,lam,mu)
        k3 = dt*self.CEOM(t+dt/2,y0+k2/2,aij,w,lam,mu)
        k4 = dt*self.CEOM(t+dt,y0+k3,aij,w,lam,mu)
        y = y0+(1/6)*(k1+2*(k2+k3)+k4)
        return y
#The Runge–Kutta–Fehlberg method
    def Integ_CEOM2(self,y0,dt,t,aij,w,lam,mu):
        k1 = dt*self.CEOM(t,y0,aij,w,lam,mu)
        k2 = dt*self.CEOM(t+dt/4,y0+k1/4,aij,w,lam,mu)
        k3 = dt*self.CEOM(t+(dt*3/8),y0+(k1*3/32)+(k2*9/32),aij,w,lam,mu)
        k4 = dt*self.CEOM(t+(dt*12/13),y0+(k1*1932/2197)-(k2*7200/2197)+(k3*7296/2197),aij,w,lam,mu)
        k5 = dt*self.CEOM(t+dt,y0+(k1*439/216)-(k2*8)+(k3*3680/513)-(k4*845/4104),aij,w,lam,mu)
        k6 = dt*self.CEOM(t+dt/2,y0-(k1*8/27)+(k2*2)-(k3*3544/2565)+(k4*1859/4104)-(k5*11/40),aij,w,lam,mu)
        z = y0+k1*16/135+k3*6656/12825+k4*28561/56430-k5*9/50+k6*2/55
        return  z

    def Integ_CEOM3(self,y0,dt,t,aij,w,lam,mu):
        k1 = dt*self.CEOM(t,y0,aij,w,lam,mu)
        k2 = dt*self.CEOM(t+dt*4/27,y0+k1*4/27,aij,w,lam,mu)
        k3 = dt*self.CEOM(t+dt*2/9,y0+(k1+3*k2)/18,aij,w,lam,mu)  
        k4 = dt*self.CEOM(t+dt/3,y0+(k1+3*k3)/12,aij,w,lam,mu)   
        k5 = dt*self.CEOM(t+dt/2,y0+(k1+3*k4)/8,aij,w,lam,mu)
        k6 = dt*self.CEOM(t+dt*2/3,y0+(13*k1-27*k3+42*k4+8*k5)/54,aij,w,lam,mu)  
        k7 = dt*self.CEOM(t+dt/6,y0+(389*k1-54*k3+966*k4-824*k5+243*k6)/4320,aij,w,lam,mu) 
        k8 = dt*self.CEOM(t+dt,y0+(-234*k1+81*k3-1164*k4+656*k5-122*k6+800*k7)/20,aij,w,lam,mu)
        k9 = dt*self.CEOM(t+dt*5/6,y0+(-127*k1+18*k3-678*k4+456*k5-9*k6+576*k7+4*k8)/288,aij,w,lam,mu)
        k10 = dt*self.CEOM(t+dt,y0+(1481*k1-81*k3+7104*k4-3376*k5+72*k6-5040*k7-60*k8+720*k9)/820,aij,w,lam,mu)
        y = y0+(41*k1+27*k4+272*k5+27*k6+216*k7+216*k9+41*k10)/840
        return y
##Sampling
    def sampling(self,w,nsamp):
        sample_x = []
        sample_v = []
        for ww in w:
            p = np.linspace(-3*np.sqrt(ww),3*np.sqrt(ww),20000)
            x = np.linspace(-3*np.sqrt(1/ww),3*np.sqrt(1/ww),20000)
            p_dis = [np.repeat(pp,1+int(1000*self.gauss(pp,ww))) for pp in p]
            x_dis = [np.repeat(xx,1+int(1000*self.gauss(xx,1/ww))) for xx in x]
            x_dis = np.asarray(self.flatt(x_dis))
            p_dis =  np.asarray(self.flatt(p_dis))
            sample_x.append(np.random.choice(x_dis,nsamp))
            sample_v.append(np.random.choice(p_dis,nsamp))
        return sample_x,sample_v

### distribution
    def gauss(self,x,d):
        return np.exp(-x**2/d)/np.pi
####flatten list
    def flatt(self,x):
        y=[]
        for l in x:
            for k in l:
                y.append(k)
        return y

