##### Generating functions
#
def MixtMod(n=2,t=1024,sigma1=1,p1=0.2,ptot=0.1):

    import numpy as np

    S = np.zeros((n,t))
    Sc_part = np.zeros((n,t))

    num_cor = np.floor(p1*ptot*t)
    num_cor = int(num_cor)
    ind_cor = randperm(t)
    ind_cor = ind_cor[0:num_cor]
    ind = np.ones((1,t))
    ind[0,ind_cor] = 0     
    ind_nocor = np.where(ind[0,:] == 1)[0] 
    
    rest = t - num_cor
    
    for r in range(0,n):
        p_active = np.floor(t*ptot*(1.-p1))
        temp = np.random.randn(1,rest)
        ind = randperm(rest)
        temp[0,ind[int(p_active+1):rest]] = 0
        cor_val = sigma1*np.random.randn(1,num_cor)
        S[r,ind_cor] = cor_val
        S[r,ind_nocor] = temp
        Sc_part[r,ind_cor] = cor_val
        
    return S,Sc_part,ind_cor,ind_nocor
    
#
def Make_Experiment_Coherent(n_s=2,n_obs=2,t_samp=1024,w=15,noise_level=40,dynamic=10,sigma1=1,p1=0.1,ptot=0.1):

    import numpy as np

    S,Sc_part,ind_cor,ind_nocor = MixtMod(n_s,t_samp,sigma1,p1,ptot)

    x = np.linspace(1,t_samp,t_samp)-t_samp/2

    kern = np.exp(-abs(x)/(w/np.log(2)))
    kern = kern/np.max(kern)

    for r in range(0,n_s):
        S[r,:] = np.convolve(S[r,:],kern,mode='same')

    val = np.power(10,(-np.linspace(1,n_s-1,n_s)/(n_s-1)*dynamic))

    A0 = np.random.randn(n_obs,n_s)
    A0 = np.dot(A0,np.diag(1./np.sqrt(np.sum(A0*A0,axis=0))))
    
    S0= np.dot(np.diag(1./np.sqrt(np.sum(S*S,axis=1))),S)
    S0 = np.dot(np.diag(val),S0)
    X0 = np.dot(A0,S0)

    N = np.random.randn(n_obs,t_samp)
    sigma_noise = np.power(10,(-noise_level/20))*np.linalg.norm(X0,ord='fro')/np.linalg.norm(N,ord='fro')
    N = sigma_noise*N

    X = X0 + N

    return X,X0,A0,S0,N,sigma_noise,kern

def randperm(n=1):

    import numpy as np

    X = np.random.randn(n)
    I = X.argsort()
    
    return I

def Gen_BG_Sources(n=2,t=1024,p=0.1):
    
    import numpy as np
    
    S = np.zeros((n,t))
    
    K = np.floor(p*t)
    
    for r in range(0,n):
        I = randperm(t)
        S[r,I[0:K]] = np.random.randn(K)
        
    return S
    
def Gen_BG_Mixtures(m=2,n=2,t=1024,p=0.1):
    
    import numpy as np
    
    S = Gen_BG_Sources(n,t,p)
    
    A = np.random.randn(m,n)
    
    A = np.dot(A,np.diag(1./np.sqrt(np.sum(A*A,axis=0))))
    
    X = np.dot(A,S)
    
    return X,A,S
    
def Gen_BG_Mixtures_Outliers(m=2,n=2,t=1024,p=0.1,po=0.01):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    S = Gen_BG_Sources(n,t,p)
    
    O = np.random.randn(m,t)
    K = np.floor(po*t*m) 
    I = np.sort(abs(np.reshape(O,(1,m*t))))
    thrd = I[0,m*t - K]
    O[abs(O) < thrd] = 0    
    
    A = np.random.randn(m,n)
    
    A = np.dot(A,np.diag(1./np.sqrt(np.sum(A*A,axis=0))))
    
    X = np.dot(A,S) + O
    
    return X,A,S,O
    
def CorrectPerm(A0,S0,A,S):

    import numpy as np
    import scipy.linalg as lng 
        
    Diff = np.dot(lng.inv(np.dot(A0.T,A0)),np.dot(A0.T,A))
    
    z = np.shape(A)
    
    for ns in range(0,z[1]):
        Diff[ns,:] = abs(Diff[ns,:])/max(abs(Diff[ns,:]))
        
    Q = np.ones(z)
    Sq = np.ones(np.shape(S))
    
    for ns in range(0,z[1]):
        Q[:,np.nanargmax(Diff[ns,:])] = A[:,ns]
        Sq[np.nanargmax(Diff[ns,:]),:] = S[ns,:]
            
    S = Sq
    A = Q
    
    return A,S
 
 ################# CODE TO COMPUTE THE MIXING MATRIX CRITERION (AND SOLVES THE PERMUTATION INDETERMINACY)    
    
def EvalCriterion(A0,S0,A,S):

    import numpy as np
    import scipy.linalg as lng 
    import copy as cp
        
    c_A0 = cp.copy(A0)
    c_A = cp.copy(A)
  
    nX = np.shape(c_A);  
    
    for r in range(0,nX[1]):
        c_A[:,r] = c_A[:,r]/(1e-24+lng.norm(c_A[:,r]))
        c_A0[:,r] = c_A0[:,r]/(1e-24+lng.norm(c_A0[:,r]))
  
    Diff = np.dot(np.linalg.inv(np.dot(c_A.T,c_A)),np.dot(c_A.T,c_A0))
                
    z = np.shape(c_A)
    
    for ns in range(0,z[1]):
        Diff[ns,:] = abs(Diff[ns,:])/max(abs(Diff[ns,:]))
        
    Q = np.ones(z)
    Sq = np.ones(np.shape(S))
    
    for ns in range(0,z[1]):
        Q[:,np.nanargmax(Diff[ns,:])] = A[:,ns]
        Sq[np.nanargmax(Diff[ns,:]),:] = S[ns,:]
            
    Diff = np.dot(np.linalg.inv(np.dot(A0.T,A0)),np.dot(A0.T,Q))
    
    for ns in range(0,z[1]):
        Diff[ns,:] = abs(Diff[ns,:])/max(abs(Diff[ns,:]))


    p = (np.sum(Diff) - z[1])/(z[1]*(z[1]-1))
    
    return p