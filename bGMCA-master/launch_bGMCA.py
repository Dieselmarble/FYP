#  Version
#    v1 - September,28 2017 - C.Kervazo - CEA Saclay 
#
#


#%%

import BSS_Utils as bu
import bGMCA_block as cgb
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import copy as cp

#%%
def launch_GMCA(n_sources=[50], blocksize=[1,2,3,4,5,7,10,13,16,18,20,28,36,45,50], numExp=1,N_MC=1,sauv=2,init=1,dataType=1,optBlock=0,folderInitSave='Results/', t_samples=1000,ptot=0.1,cd=1,SNR=120,verPos=0,J=0):
    '''Usage:
    launch_GMCA(numExp=1,N_MC=1,sauv=2,init=1,optBlock=0,folderInitSave=0, t_samples=0,ptot=0.1,cd=1,colSize=0,palm=0,SNR=120,optWhite=0,verPos=0,J=0)
    
    Inputs:
    n_sources : list of intergers (list containing the number of sources of the experiments)
    blocksize : list of intergers (list of the different block sizes to be tested. For a given number of sources n in the n_sources list, valid if the blocksizes are in [1,n])
    numExp : int (number of the experiment (used for the seed of the random data generation))
    N_MC : int (number of Monte-Carlo simulations)
    sauv : int (0 : no saving, 1 : save the results, 2 : save all including data
    init : boolean (0: load some already existing X,X0,A0,S0,N data matrices, 1 : creation of new random data)
    dataType : int in {1,2,3} (1: Bernouilli Gaussian data, 2: Generalized Gaussian data (approximately sparse), 3: realistic 1H NMR data)
    optBlock : int in {1,2,3} (0 : blocks are crated randomly, 1 : creation according to the "correlation" of the sources, 2 : according to the angle between the columns of A)
    folderInitSave : string (name of the folder use to save or load the data. This folder must contain 2 subfolders 'Data' and 'MatricesAS')
    t_samples : int (number of samples for each source)
    ptot : float in (0,1) (sparsity of the sources. More sparse for low values:
                            if dataType = 1: ptot = number of non-zero coefficients
                            if dataType = 2: ptot = parameter of Generalized Gaussian: 1 => Laplacian, 2 = Gaussian
                            if dataType = 3: ptot is not used)
    cd : int >= 1 (condition number of A)
    SNR : float (signal to noise ratio. Do not use int)
    verPos : int in {0,1,2,3} (0 : no non-negativity, verPos = 1 : wavelets + non-negativity in the direct domain, verPos = 2 : wavelets without non-negativity, verPos = 3 : non-negativity)
    J : int (number of wavelets scales if applicable)
    
    Outputs:
    C_hyb : len(blocksize) x len(n_sources) x N_MC array (contains the C_A values for all the blocksizes and different number of sources for each monte carlo experiment)
    S_gmca : n x t array (found sources. Relevant only if len(blocksize) = len(n_sources) = N_MC = 1)
    A_gmca : m x n array (found mixing matrix. Relevant only if len(blocksize) = len(n_sources) = N_MC = 1)'''

    nitGMCA = 1000 # Number of iterations, should be fixed so that nit > 100*n_sources/blocksize
    nitPALM = nitGMCA*10
    



    pltRes = 1 # Boolean value




    # bGMCA is randomly initialized
    initGMCA = 1
    

    
    if init == 0 and sauv == 2:
        print('Writing and reading of data at the same time')
        return
    


    # Creation of strings for the file name    
    cdStr = cd


    if ptot == 1:
        ptotStr = '1'
    else:
        ptotStr = '0_' + str(ptot-int(ptot))[2:]

           
    sourcesStr = str(n_sources)
    
    # Initializations
    C_hyb = np.zeros((len(blocksize),len(n_sources),N_MC))
    boolVal = 1
    
    
    
    for It_MC in range(0,N_MC): # Loop for performing N_MC monte carlo simulations
        valSeed = numExp*N_MC + It_MC + 1
        np.random.seed(valSeed)
        print('Value of the seed %s'%(valSeed))
        for R_t in range(0,len(n_sources)): # Loop on the different number of sources to be tested
            
            n = n_sources[R_t]
    
            # Generating some random sparse sources with random mixing matrices
            # in the determined case (m=n)
            if init ==1:
                if dataType == 1:
                    X,X0,A0,S0,N = bu.Make_Experiment_Exact(n_s=n,n_obs=n,t_samp=t_samples,noise_level=SNR,dynamic=0,ptot=ptot,cd=cd)       
                elif dataType == 2:
                    X,X0,A0,S0,N = bu.Make_Experiment_GeneralizedGaussian(n_s=n,n_obs=n,t_samp=t_samples,noise_level=120,dynamic=0,ptot=1,cd=cd,alpha=ptot)       
                elif dataType == 3:
                    X,X0,A0,S0,N = bu.Make_Experiment_HMR(n_s=n,n_obs=8*n,t_samp=t_samples,noise_level=SNR,standDev = 3,peak_width=3)


            if init == 0: # If we do not create new data matrices but rather use some already saved matrices
                X = sio.loadmat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_X_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,0,1,SNR,3,J),{'donnees':X})
                X = X['donnees']
                X0 = sio.loadmat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_X0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,0,1,SNR,3,J),{'donnees':X0})
                X0 = X0['donnees']
                A0 = sio.loadmat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_A0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,0,1,SNR,3,J),{'donnees':A0})
                A0 = A0['donnees']
                S0 = sio.loadmat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_S0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,0,1,SNR,3,J),{'donnees':S0})
                S0 = S0['donnees']
                N = sio.loadmat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_N_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,0,1,SNR,3,J),{'donnees':N})                
                N = N['donnees']
            

                    
            if sauv == 2 and init ==  1: # The data matrices are saved 
                sio.savemat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_X_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':X})
                sio.savemat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_X0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':X0})
                sio.savemat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_A0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':A0})
                sio.savemat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_S0_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':S0})
                sio.savemat(folderInitSave + 'Data/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_N_exp%s_MC%s_sourceN%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':N})
            
            for R_p in range(0,len(blocksize)): # Loop on the different block sizes to be tested
                minibatch = blocksize[R_p] 
                
                if minibatch > n:
                    C_hyb[R_p,R_t,It_MC] = -1
                    
                else:
                    if verPos == 0:# No wavelet and no non-negativity constraint
                        g_S_gmca,g_A_gmca,exception = cgb.bGMCA(X,n,Init=initGMCA,mints=1,nmax=nitGMCA,L0=0,verb=1,blocksize=minibatch,optBlock=optBlock,J=0)
                        
                    elif verPos == 1:# Wavelet and non-negativity in the direct domain
                        g_S_gmca,g_A_gmca,exception = cgb.bGMCA_NMF_ondt_naif(X,n,Init=initGMCA,mints=1,nmax=nitGMCA,L0=0,verb=1,blocksize=minibatch,optBlock=optBlock,J=J)                       

                    elif verPos == 2:# Wavelets without non-negativity
                        g_S_gmca,g_A_gmca,exception,Sw = cgb.bGMCA(X,n,Init=initGMCA,mints=1,nmax=nitGMCA,L0=0,verb=1,blocksize=minibatch,optBlock=optBlock,J=J)

                    elif verPos == 3: # Non-negativity in the direct domain
                        g_S_gmca,g_A_gmca,exception = cgb.bGMCA(X,n,Init=initGMCA,mints=1,nmax=nitGMCA,L0=0,verb=1,blocksize=minibatch,optBlock=optBlock,J=0,optPos=1)
                        
                        
                    print('End of GMCA stage')                        
                    S_gmca = g_S_gmca
                    A_gmca = g_A_gmca
                        
                    
                    
                    if exception != 0:
                        boolVal = 0
                        print('Exception = %s' %exception)
                                    
                                    
                                    
                                    
                    if verPos == 0:# No wavelet and no non-negativity constraint
                        g_S_gmca,g_A_gmca = cgb.PALM_NMF_MainBlock(X,n,A=cp.deepcopy(A_gmca),S=cp.deepcopy(S_gmca),kend=1,nmax=nitPALM,L0=0,blocksize=minibatch,tol=1e-12)  
                    elif verPos == 1:# Wavelet and non-negativity in the direct domain
                        g_S_palm,g_A_gmca = cgb.PALM_NMF_MainBlock_prox(X,n,A=cp.deepcopy(A_gmca),S=cp.deepcopy(S_gmca),kend=1,nmax=nitPALM,L0=0,verb=0,blocksize=minibatch,tol=1e-12,J=J)  
                    elif verPos == 2:# Wavelets without non-negativity
                        g_S_gmca,g_A_gmca = cgb.PALM_NMF_MainBlock(X,n,A=cp.deepcopy(A_gmca),S=cp.deepcopy(Sw),kend=1,nmax=nitPALM,L0=0,blocksize=minibatch,tol=1e-12,J=J)  
                    elif verPos == 3: # Non-negativity in the direct domain
                        g_S_gmca,g_A_gmca = cgb.PALM_NMF_MainBlock(X,n,A=cp.deepcopy(A_gmca),S=cp.deepcopy(S_gmca),kend=1,nmax=nitPALM,L0=0,blocksize=minibatch,tol=1e-12,optPos=1)

                    
                    S_hyb = g_S_gmca
                    A_hyb = g_A_gmca
                    
                    
                    
                    try:
                        C_hyb[R_p,R_t,It_MC] = bu.EvalCriterion(A0,S0,cp.deepcopy(A_hyb),cp.deepcopy(S_hyb))
                        
                    except ValueError:
                        C_hyb[R_p,R_t,It_MC] = np.log10(-1)
                            
                    
                    
                        
                    
                    if sauv > 0:# We save the A and S matrices
                        sio.savemat(folderInitSave + 'MatricesAS/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_S_gmca_exp%s_MC%s_sourceN%s_blocksizeNumber%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,R_p,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':S_gmca})
                        sio.savemat(folderInitSave + 'MatricesAS/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_A_gmca_exp%s_MC%s_sourceN%s_blocksizeNumber%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,R_p,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':A_gmca})
                                          
                        sio.savemat(folderInitSave + 'MatricesAS/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_S_hyb_exp%s_MC%s_sourceN%s_blocksizeNumber%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,R_p,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':S_hyb})
                        sio.savemat(folderInitSave + 'MatricesAS/GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_A_hyb_exp%s_MC%s_sourceN%s_blocksizeNumber%s_optBlock%s_init%s_SNR%s_optWhite%s_verPos%s_J%s'%(sourcesStr,numExp,It_MC,R_t,R_p,optBlock,initGMCA,SNR,0,verPos,J),{'donnees':A_hyb})
                                                              
        
    if sauv > 0:# The block sizes and the corresponding C_A values are saved
        sio.savemat(folderInitSave + 'GMCA_cluster_GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_optBlock%s_blocksize_init%s_SNR%s_optWhite%s_verPos%s_J%s_exp%s'%(sourcesStr,optBlock,initGMCA,SNR,0,verPos,J,numExp),{'donnees':blocksize})
              
        sio.savemat(folderInitSave + 'GMCA_cluster_GMCA_cond%s_%sech_%sit_ptot'%(cdStr,t_samples,nitGMCA) + ptotStr + '_sources%s_optBlock%s_C_hyb_init%s_SNR%s_optWhite%s_verPos%s_J%s_exp%s'%(sourcesStr,optBlock,initGMCA,SNR,0,verPos,J,numExp),{'donnees':C_hyb})

                           
    if pltRes > 0:
        font = {'weight' : 'bold', 'size'   : 20}
                   
        plt.figure(1)
        plt.hold(True)
        plt.plot(blocksize,-10*np.log10(np.nanmedian(C_hyb,axis=2)),linewidth = 4)
        
        
        plt.figure(1)
        plt.xlabel('Block sizes',**font)            
        plt.ylabel('Mixing matrix criterion in log scale',**font)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.title('Mixing matrix criterion using different block sizes')
       

        
        if boolVal == 0:
            print('WARNING, SVD DID NOT CONVERGE FOR SOME EXPERIMENTS')
        else:
            print('Correct termination of the algorithm')

    
        
    return C_hyb,S_gmca,A_gmca