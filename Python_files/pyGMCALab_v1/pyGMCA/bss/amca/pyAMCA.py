#
#  Usage:
#    S,A,PiA,WeiA,wL1 = pam.AMCA(X0,n,mints=3,nmax=250,q_f = 0.1,AMCA=1,rL1=1)
#
#
#  Inputs:
#   X0    : m x t array (input data, each row corresponds to a single observation)
#   n     : scalar (number of sources to be estimated)
#   mints : scalar (final value of the k-mad thresholding)
#   nmax  : scalar (number of iterations)
#   q_f   : scalar (final value of the weighting parameter) 
#   AMCA  : boolean (if set to 1, runs AMCA; otherwise GMCA is performed)   
#   rL1   : boolean (if set to 1, applies reweighting L1)
#   UseP  : boolean (if set to 1, thresholds as selected based on the selected source support: a fixed number of active entries is selected at each iteration)
#   Init  : scalar (if set to 0: PCA-based initialization, if set to 1: random initialization)
#   verb   : boolean (if set to 1, in verbose mode)
#
#  Outputs:
#   S    : n x t array (estimated sources)
#   A    : m x n array (estimated mixing matrix)
#
#  Description:
#    Computes the sources and the mixing matrix with AMCA or GMCA
#
#  Example:
#     S,A = pam.AMCA(X,n,mints=0,nmax=500,AMCA=0,UseP=1) will perform GMCA assuming that the data are noiseless
#
#  Version
#    v1 - April,14 2015 - J.Bobin - CEA Saclay 
#
#
################# DEFINES THE MEDIAN ABSOLUTE DEVIATION    

def mad(xin = 0):

    import numpy as np

    z = np.median(abs(xin - np.median(xin)))/0.6735
    
    return z

################# AMCA Main function    

def AMCA(X,n,maxts = 0,mints=3,nmax=100,q_f = 0.1,AMCA=0,rL1=0,L0=0,UseP=1,verb=0,Init=0):

    import numpy as np
    import scipy.linalg as lng 
    import copy as cp
    import matplotlib.pyplot as plt
#
    nX = np.shape(X);
    m = nX[0];t = nX[1]    
    Xw = cp.copy(X)
#
#---- Main loop
#Initialization could largely be improved -> it should be data-dependent
    if verb:
        print("Initializing ...")
    if Init == 0:
        R = np.dot(Xw,Xw.T)
        D,V = lng.eig(R)
        A = V[:,0:n]
    if Init == 1:
        A = np.random.randn(m,n)    
    for r in range(0,n):
        A[:,r] = A[:,r]/lng.norm(A[:,r]) # - We should do better than that
#   
    S = np.dot(A.T,Xw);   
    maxts = 0
    for r in range(0,n):
        indNZ = np.where(abs(S[r,:]) > 0)[0]
        gg = np.max(abs(S[r,indNZ]))/mad(S[r,indNZ])
        maxts = np.max([gg,maxts]) 
    maxts = np.max([maxts,10.])        
    S,A,PiA,WeiA,wL1 = GMCA_MainBlock(Xw,n,A,S,q_f=q_f,kend = mints,maxts=maxts,nmax=nmax,AMCA=AMCA,rL1=rL1,L0=L0,UseP=UseP,verb=verb);
#    
    return S,A

################# AMCA internal code (Weighting the sources only)

def GMCA_MainBlock(X=0,n=0,A=0,S=0,q_f=0.1,kend=3,maxts=5,nmax=100,AMCA=0,rL1=0,L0=1,UseP=0,verb=0):
#--- Import useful modules 
    import numpy as np
    import scipy.linalg as lng 
    import copy as cp
    import scipy.io as sio
    import time
#--- Init
    n_X = np.shape(X)    
    n_S = np.shape(S)
    W = np.ones((1,n_S[1]));    
    wL1 = np.ones(np.shape(S));    
    k = maxts
    dk = (k-kend)/(nmax-1);
    alpha = 1;
    dalpha = (alpha-q_f)/nmax;
    perc = 1./nmax
    Aold = cp.deepcopy(A)
#    
    Go_On = 1
    it = 1
#   
    if verb:
        print("Starting main loop ...")    
        print(" ")
        print("  - Final k: ",kend)
        print("  - Maximum number of iterations: ",nmax)
        if UseP:
            print("  - Using support-based threshold estimation")
        if rL1:
            print("  - Using reweighting L1")
        if L0:
            print("  - Using L0 norm rather than L1")
        if AMCA:
            print("  - AMCA rather than GMCA")
        print(" ")
        print(" ... processing ...")
        start_time = time.time()
#--- Main loop
    while Go_On:
        it += 1
        if it == nmax:
            Go_On = 0
    #for it in range(0,nmax):
        #--- Estimate the sources
        sigA = np.sum(A*A,axis=0)
        indS = np.where(sigA > 0)[0]
        if np.size(indS) > 0: 
            Ra = np.dot(A[:,indS].T,A[:,indS])
            Ua,Sa,Va = np.linalg.svd(Ra)
            cd_Ra = np.min(Sa)/np.max(Sa)
            if cd_Ra > 1e-6:
                iRa = np.dot(Va.T,np.dot(np.diag(1./Sa),Ua.T))
                piA = np.dot(iRa,A[:,indS].T)    
                S[indS,:] = np.dot(piA,X)
            if cd_Ra < 1e-6:
                #--- iterative update
                La = np.max(Sa)
                for it_A in range(0,10):
                    S[indS,:] = S[indS,:] + 1/La*np.dot(A[:,indS].T,X - np.dot(A[:,indS],S[indS,:]))
                    # We could threshold as well
            Stemp = S[indS,:]
            Sref = cp.copy(S)
            for r in range(np.size(indS)):
                St = Stemp[r,:]
                indNZ = np.where(abs(St) > kend*mad(St))[0]
                thrd = mad(St[indNZ])
                if UseP == 0:
                    thrd = k*thrd
                if UseP == 1:
                    Kval = np.min([np.floor(perc*it*len(indNZ)),n_S[1]-1.])
                    I = abs(St[indNZ]).argsort()[::-1]
                    Kval = np.min([np.max([Kval,n]),len(I)-1.])
                    thrd = abs(St[indNZ[I[int(Kval)]]])
                St[(abs(St)/wL1[indS[r],:] < thrd)] = 0 #--- Be careful
                indNZ = np.where(abs(St) > thrd)[0] 
                if L0 == 0:
                    St[indNZ] = St[indNZ] - thrd*wL1[indS[r],indNZ]*np.sign(St[indNZ])
                if rL1 > 0:
                    wL1[indS[r],:] = np.exp(-abs(St)/(np.median(abs(St[indNZ])) + mad(abs(St[indNZ]))))+1e-12 #---- Should be improve
                Stemp[r,:] = St    
            S[indS,:] = Stemp
            k = k - dk
            # --- Updating the weighting
            Ns = np.sqrt(np.sum(S*S,axis=1))
            IndS = np.where(Ns > 0)[0]
            if AMCA > 0:   
                if len(IndS) > 0:
                    
                    Sref[IndS,:] = np.dot(np.diag(1./Ns[IndS]),Sref[IndS,:])
                    W = np.power(np.sum(np.power(abs(Sref[indS,:]),alpha),axis=0),1./alpha)
                    ind = np.where(W > 0)[0]
                    jind = np.where(W == 0)[0]
                    W[ind] = 1./W[ind];        
                    if len(jind) > 0:
                        W[jind] = np.median(W[ind]) 
                    alpha = alpha-dalpha
                W /= np.max(W)
        # --- Updating the mixing matrix
        Ns = np.sqrt(np.sum(S*S,axis=1))
        indA = np.where(Ns > 0)[0]
        if len(indA) > 0:
            Sr = cp.deepcopy(S)
            for r in range(0,n):
                Sr[r,:] = Sr[r,:]*W
            Rs = np.dot(S[indA,:],Sr[indA,:].T)
            Us,Ss,Vs = np.linalg.svd(Rs)
            cd_Rs = np.min(Ss)/np.max(Ss)
            if cd_Rs > 1e-6:
                piS = np.dot(Sr[indA,:].T,np.linalg.inv(Rs));
                A[:,indA] = np.dot(X,piS)
                A = np.dot(A,np.diag(1./(1e-24 + np.sqrt(np.sum(A*A,axis=0)))));
            if cd_Rs < 1e-6:
                #--- iterative update
                Ls = np.max(Ss)
                for it_A in range(0,10):
                    A[:,indA] = A[:,indA] + 1/Ls*np.dot(X - np.dot(A[:,indA],S[indA,:]),Sr[indA,:].T)
                    A[:,indA] = np.dot(A[:,indA],np.diag(1./(1e-24 + np.sqrt(np.sum(A[:,indA]*A[:,indA],axis=0)))));
        DeltaA =  np.linalg.norm(Aold - A)
        if DeltaA < 1e-10:
            if it > 100:
                Go_On = 0
        Aold = cp.deepcopy(A)
    if verb:
        elapsed_time = time.time() - start_time
        print("Stopped after ",it," iterations, in ",elapsed_time," seconds")
    #
    return S,A,piA,W,wL1
#
################# END OF AMCA