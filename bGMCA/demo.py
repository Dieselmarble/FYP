# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:24:53 2017

@author: ckervazo
"""

import launch_bGMCA as lgc
#%% Simple experience
N_MC = 1# only one Monte Carlo experiment is performed
numExp = 2# used for setting the random seed for the data generation
folderInitSave='ResultsDedale/' # folder to write the results. It must contain 2 subfolders « data » and « matricesAS ».
optBatch = 0# correspond to random choices for the blocks
ptot = 0.1# proportion of non-zero coefficients
cd = 1# condition number of the mixing matrix

dataType = 1 #The data is created a Bernouilli Gaussian distribution: a proportion of ptot coefficients is taken non-zeros and have a Gaussian amplitude.
verPos = 0 #Corresponds to a simple sparsity constraint in the direct domain

blocksize = [1,2,3,4,5,7,10,13,16,18,20] # Block sizes to be tested
n_sources = [20] #Number of sources. 20 sources is already a quite high number of sources. Take 50+ sources for really large-scale examples.

C_gmca,S,A = lgc.launch_GMCA(n_sources=n_sources,blocksize=blocksize,numExp=numExp,N_MC=N_MC,dataType=dataType,sauv=2,optBlock=optBatch,folderInitSave=folderInitSave,ptot=ptot,verPos=verPos,cd=cd)

#%% Experience with wavelets and non negativity
#N_MC = 1# only one Monte Carlo experiment is performed
#numExp = 2# used for setting the random seed for the data generation
#folderInitSave='ResultsDedale/' # folder to write the results. It must contain 2 subfolders « data » and « matricesAS ».
#optBatch = 0# correspond to random choices for the blocks
#
#dataType = 3 #The data comes from saved realistic sources. All the matrices are non-negative.
#verPos = 1 #Enforces the non-negativity and sparsity constraint.
#J=2 # Number of wavelet scales.
#
#blocksize = [1,2,3,4,5,7,10,13,16,18,20] # Block sizes to be tested
#n_sources = [20] #Number of sources. 20 sources is already a quite high number of sources. Take 50+ sources for really large-scale examples.
#
#C_gmca,S,A = lgc.launch_GMCA(n_sources=n_sources,blocksize=blocksize,numExp=numExp,N_MC=N_MC,dataType=dataType,sauv=2,optBlock=optBatch,folderInitSave=folderInitSave,verPos=verPos,J=J)

#%% Experience about the impact of the condition number
#cdTab = [2,6,20]
#N_MC = 1
#
#numExp = 2# used for setting the random seed for the data generation
#optBatch = 0# correspond to random choices for the blocks
#
#dataType = 1 #The data is created a Bernouilli Gaussian distribution: a proportion of ptot coefficients is taken non-zeros and have a Gaussian amplitude.
#verPos = 0 #Corresponds to a simple sparsity constraint in the direct domain
#
#blocksize = [1,2,3,4,5,7,10,13,16,18,20] # Block sizes to be tested
#n_sources = [20] #Number of sources. 20 sources is already a quite high number of sources. Take 50+ sources for really large-scale examples.
#
#for ii in range(len(cdTab)):
#    folderInitSave = 'ResultsDedale/%scd/'%(cdTab[ii])
#    C_gmca,S,A = lgc.launch_GMCA(n_sources=n_sources,blocksize=blocksize,numExp=numExp,N_MC=N_MC,dataType=dataType,sauv=2,optBlock=optBatch,folderInitSave=folderInitSave,verPos=verPos,cd=cdTab[ii])
    
#%% Experience about the impact of the proportion of non zero coefficients
#ptotTab = [0.1,0.2,0.3]
#N_MC = 1
#
#numExp = 2# used for setting the random seed for the data generation
#optBatch = 0# correspond to random choices for the blocks
#
#dataType = 1 #The data is created a Bernouilli Gaussian distribution: a proportion of ptot coefficients is taken non-zeros and have a Gaussian amplitude.
#verPos = 0 #Corresponds to a simple sparsity constraint in the direct domain
#
#blocksize = [1,2,3,4,5,7,10,13,16,18,20] # Block sizes to be tested
#n_sources = [20] #Number of sources. 20 sources is already a quite high number of sources. Take 50+ sources for really large-scale examples.
#
#for ii in range(len(ptotTab)):
#    ptotStr = '0_' + str(ptotTab[ii]-int(ptotTab[ii]))[2:]
#    folderInitSave = 'ResultsDedale/ptot' + ptotStr + '/'
#    C_gmca,S,A = lgc.launch_GMCA(n_sources=n_sources,blocksize=blocksize,numExp=numExp,N_MC=N_MC,dataType=dataType,sauv=2,optBlock=optBatch,folderInitSave=folderInitSave,verPos=verPos,ptot=ptotTab[ii])