"""
This code provides some examples of the AMCA algorithms

@author: J.Bobin
@date: April, 14 2015

"""
import sys
sys.path.insert(0, '/Users/kevin/Documents/MATLAB/Python_files/pyGMCALab_v1')
from  pyGMCA.bss.amca import BSS_Utils as bu
from  pyGMCA.bss.amca import pyAMCA as pam
import matplotlib.pyplot as plt
import numpy as np
#
N_MC = 25 #- Number of Monte-Carlo simulations
#
print("Comparing GMCA and AMCA when sparsity is enforced in the sample domain")
#%%
print("Evolution of the mixing matrix criterion when the sparsity level varies")
#
n = 4
m = 4
t = 2048
p1 = 0.1
sigma1 = 4
pval = np.power(10,np.linspace(-2,-0.7,5))
#
C_G_sparse_level = np.zeros((5,N_MC))
C_A_sparse_level = np.zeros((5,N_MC))
#
for It_MC in range(0,N_MC):
    for R_p in range(0,5):
        p = pval[R_p]
        X,X0,A0,S0,N,sigma_noise,kern = bu.Make_Experiment_Coherent(n_s=n,n_obs=m,t_samp=t,w=0.1,noise_level=120,dynamic=0,sigma1=sigma1,p1=p1,ptot=p)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=0,rL1=0,L0=0,UseP=1)
        C_G_sparse_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=1,rL1=0,L0=0,UseP=1)
        C_A_sparse_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
#
#%% PLOTTING THE RESULTS
plt.figure(0)
pval = np.power(10,np.linspace(-2,-0.7,5))
plt.title('Mixing matrix criterion as a function of the sparsity level')
tempG = np.median(C_G_sparse_level,1)
tempA = np.median(C_A_sparse_level,1)
Mrange = 1.5*np.max([np.max(tempA),np.max(tempG)])
mrange = 0.5*np.min([np.min(tempA),np.min(tempG)])
plt.semilogy(pval,tempG,'kd')
plt.semilogy(pval,tempA,'ro')
plt.axis([0,0.2, 1e-10,1])
#
print("Experiment complete")
#%%
print("Evolution of the mixing matrix criterion when the correlation level varies")
#
n = 4
m = 4
t =2048
p1 = 0.1
sigma1 = 4
p = 0.05
pval = np.power(10,np.linspace(-2,np.log10(0.95),10))
#
C_G_corr_level = np.zeros((10,N_MC))
C_A_corr_level = np.zeros((10,N_MC))
#
for It_MC in range(0,N_MC):
    for R_p in range(0,10):
        p1 = pval[R_p]
        X,X0,A0,S0,N,sigma_noise,kern = bu.Make_Experiment_Coherent(n_s=n,n_obs=m,t_samp=t,w=0.1,noise_level=120,dynamic=0,sigma1=sigma1,p1=p1,ptot=p)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=0,rL1=0,L0=0,UseP=1)
        C_G_corr_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=1,rL1=0,L0=0,UseP=1)
        C_A_corr_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
#
#%% PLOTTING THE RESULTS
plt.figure(1)
pval = np.power(10,np.linspace(-2,np.log10(0.95),10))
plt.title('Mixing matrix criterion as a function of the correlation level')
tempG = np.median(C_G_corr_level,1)
tempA = np.median(C_A_corr_level,1)
Mrange = 1.5*np.max([np.max(tempA),np.max(tempG)])
mrange = 0.5*np.min([np.min(tempA),np.min(tempG)])
plt.semilogy(pval,tempG,'kd')
plt.semilogy(pval,tempA,'ro')
plt.axis([np.min(pval),np.max(pval), 1e-12,1])
#
print("Experiment complete")
#%%
print("Evolution of the mixing matrix criterion when the dynamic varies")
#
n = 4
m = 4
t = 2048
p1 = 0.1
sigma1 = 4
p = 0.05
pval = np.linspace(0,4,10)
#
C_G_dyn_level = np.zeros((10,N_MC))
C_A_dyn_level = np.zeros((10,N_MC))
#
for It_MC in range(0,N_MC):
    for R_p in range(0,10):
        dyna = pval[R_p]
        X,X0,A0,S0,N,sigma_noise,kern = bu.Make_Experiment_Coherent(n_s=n,n_obs=m,t_samp=t,w=0.1,noise_level=120,dynamic=dyna,sigma1=sigma1,p1=p1,ptot=p)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=0,rL1=0,L0=0,UseP=1)
        C_G_dyn_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
        g_S,g_A = pam.AMCA(X0,n,mints=0,nmax=500,q_f = 0.1,AMCA=1,rL1=0,L0=0,UseP=1)
        C_A_dyn_level[R_p,It_MC] = bu.EvalCriterion(A0,S0,g_A,g_S)
#
#%% PLOTTING THE RESULTS
plt.figure(2)
pval = np.linspace(0,6,10)
plt.title('Mixing matrix criterion as a function of the dynamic')
tempG = np.median(C_G_dyn_level,1)
tempA = np.median(C_A_dyn_level,1)
Mrange = 1.5*np.max([np.max(tempA),np.max(tempG)])
mrange = 0.5*np.min([np.min(tempA),np.min(tempG)])
plt.semilogy(pval,tempG,'kd')
plt.semilogy(pval,tempA,'ro')
plt.axis([np.min(pval),np.max(pval), 1e-12,1])
#
print("Experiment complete")
