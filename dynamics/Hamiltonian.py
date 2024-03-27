import numpy as np

class hamiltonian(object):
    def __init__(self,epsilon,Jcou,lamb,mu,w):
        super().__init__()
        self.epsilon =  epsilon
        self.Jcou = Jcou
        self.lamb = lamb
        self.mu = mu
        self.w = w

###electronic Hamiltonian of the type H= \sum_i \epsilon*\hat{a}^{\dagger}_{i}*\hat{a}_{i}+\sum_{i,j}Jcou_{i,j}*\hat{a}^{\dagger}_{i}*\hat{a}_{j}
    def Hel(self):
        return self.epsilon+self.Jcou#+mu
#####electronic-classic Hamiltonian considering spin-orbit H_{el-c}= H_{SO}^{el-c}+H_{el-c}
    def Helc(self,x):
        self.Hcla = np.zeros((len(self.mu),len(self.mu)))
        for i in range(len(x)):
            self.Hcla += self.w[i]*self.lamb[i]*self.mu*x[i]
        #self.Hcla = self.Hcla*(np.triu(np.tri(len(self.Jcou),len(self.Jcou),1))+
        #                       np.tri(len(self.Jcou),len(self.Jcou),-1)-np.eye(len(self.Jcou)))
        return self.Hcla
