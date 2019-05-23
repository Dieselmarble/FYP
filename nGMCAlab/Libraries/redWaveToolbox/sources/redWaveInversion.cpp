// redWaveInversion.cpp - This file is part of the redWaveToolbox.
// This software aims at performing redundant wavelet transformations.
// Copyright 2014 CEA
// Contributor : Jeremy Rapin (jeremy.rapin@cea.fr)
// Created on 14/7/2014, last modified on 16/7/2014
// 
// This software is governed by the CeCILL license under French law and
// abiding by the rules of distribution of free software. You can use,
// modify and/ or redistribute the software under the terms of the CeCILL
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
// 
// As a counterpart to the access to the source code and rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty and the software's author,  the holder of the
// economic rights,  and the successive licensors have only limited
// liability.
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean that it is complicated to manipulate,  and that also
// therefore means that it is reserved for developers and experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL license and that you accept its terms.



#define RW_MATLAB_INTERFACE
#include "redWaveTools.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //inputs:
    //- x: input data
    //- AtA
    //- AtY
    //- threshold
    //- parameter structure with optional fields:
    //  - qmf: filter (default: Daubechies-4)
    //  - dimensions: vector containing the wavelet dimension(s) (default: )
    //  - L: number of scales (default: 3)
    //  - MaximumIteration: number of iterations (default: 20)
    //  - nonNegative: find non-negative sources only (default: 0)
    //  - isometric: using isometric transform or noise normalized transform (optional, default = 0)
    
    
    ///////////////////////////////////////////////
    //     check for correct number of inputs     //
    ////////////////////////////////////////////////
    if (nrhs<4) {
        RW_PRINT("ERROR: At least 4 input parameters are required:\n- the input matrix.\n- A' * A.\n- A' * Y.\n- the thresholds.\n");
        RW_PRINT("The optional 5th input parameter structure can include fields:\n");
        RW_PRINT("   - qmf: filter (default: Daubechies-4)\n");
        RW_PRINT("   - dimensions: vector containing the wavelet dimension(s) (default: all but dimension 1)\n");
        RW_PRINT("   - L: number of scales (default: 3)\n");
        RW_PRINT("   - MaximumIteration: number of iterations (default: 24)\n");
        RW_PRINT("   - nonNegative: find non-negative sources only (default: 0)\n");
        RW_PRINT("   - isometric: using isometric transform or noise normalized transform (default = 0)\n");
        return;
    }
    bool error = false;
    
    ////////////////////////////
    //     input matrices     //
    ////////////////////////////
    
    NdimData x0 = matlabArrayToNdimData(prhs[0]);
    NdimData AtA = matlabArrayToNdimData(prhs[1]);
    NdimData AtY = matlabArrayToNdimData(prhs[2]);
    NdimData thresholds = matlabArrayToNdimData(prhs[3]);
    
    /////////////////////////////
    //     optional fields     //
    /////////////////////////////
    
    //default value of parameters fields
    int verbose = 0;
    int maximumIteration = 24;
    int isometric = 0;
    int L = 3;
    int nonNegative = 0;
    vectorData<double>* filter = NULL;
    vectorData<int>* waveDims = NULL;

    
    if (nrhs > 4) {
        const mxArray *param = prhs[4];
        
        //verbose
        int field_num = mxGetFieldNumber(param, "verbose");
        if (field_num >= 0)
            verbose = (int) *mxGetPr(mxGetFieldByNumber(param, 0, field_num));
        
        // number of scales
        field_num = mxGetFieldNumber(param, "L");
        if (field_num >= 0)
            L = (int) *mxGetPr(mxGetFieldByNumber(param, 0, field_num));
        else
            if (verbose)
                RW_PRINT("Number of scales set to default value (%i)\n", L);
        
        // number of scales
        field_num = mxGetFieldNumber(param, "isometric");
        if (field_num >= 0)
            isometric = (int) *mxGetPr(mxGetFieldByNumber(param, 0, field_num));
        else
            if (verbose)
                RW_PRINT("Isometric mode ''isometric'' set to default value (%i).\n", isometric);
        
        // maximum iteration
        field_num = mxGetFieldNumber(param, "MaximumIteration");
        if (field_num >= 0)
            maximumIteration = (int) *mxGetPr(mxGetFieldByNumber(param, 0, field_num));
        else
            if (verbose)
                RW_PRINT("Maximum number of iterations ''MaximumIteration'' set to default value (%i).\n", maximumIteration);
        
        // number of scales
        field_num = mxGetFieldNumber(param, "nonNegative");
        if (field_num >= 0)
            nonNegative = (int) *mxGetPr(mxGetFieldByNumber(param, 0, field_num));
        else
            if (verbose)
                RW_PRINT("Non-negative mode ''nonNegative'' set to default value (%i).\n", nonNegative);
        
        //filter
        field_num = mxGetFieldNumber(param, "qmf");
        if (field_num >= 0) {
            mxArray *tmp = mxGetFieldByNumber(param, 0, field_num);
            filter = new vectorData<double>(max(mxGetN(tmp), mxGetM(tmp)), mxGetPr(tmp));
        }
        else
            if (verbose)
                RW_PRINT("Filter ''qmf'' set to default value (Daubechies-4).\n");
        
        //wavelet dimensions
        field_num = mxGetFieldNumber(param, "dimensions");
        if (field_num >= 0) {
            mxArray *tmp = mxGetFieldByNumber(param, 0, field_num);
            vectorData<double> waveDimsIn(max(mxGetN(tmp), mxGetM(tmp)), mxGetPr(tmp));
            waveDims = new vectorData<int>(max(mxGetN(tmp), mxGetM(tmp)));
            for (int k = 0; k < waveDimsIn.getLength(); ++k)
                waveDims->at(k) = (int) waveDimsIn.at(k) - 1;
        }
        else
            if (verbose)
                RW_PRINT("Wavelet dimensions ''dimensions'' set to default value (all but 1).\n");
        
        
    }
    

    
    if (!waveDims) {
        waveDims = new vectorData<int>(x0.getNdim() - 1);
        for (int k = 0; k < x0.getNdim() - 1; ++k)
            waveDims->at(k) = k + 1;
    }
    if (!filter) {
        filter = new vectorData<double>(4);
        filter->at(0) = 0.482962913145;
        filter->at(1) = 0.836516303738;
        filter->at(2) = 0.224143868042;
        filter->at(3) = -0.129409522551;
    }
    
    
    
    //////////////
    // Wavelets //
    //////////////
    
    // check number of wavelet levels
    if (L < 0) {
        RW_PRINT("ERROR: The number of levels, L, must be a non-negative integer.\n");
        error = true;
    }
    
    //create wavelets
    RedWave W(*waveDims, *filter, L, isometric);
    // do not forget to initialize the sizes before using the wavelets
    if (W.setSizes(x0, true, true))
        error = true;
    
    
    if (!W.isWaveSize(thresholds)) {
        error = true;
        RW_PRINT("ERROR: thresholds matrix must be of wavelet size.\n");
    }
    
    
    
    ////////////////////
    //     output     //
    ////////////////////
    plhs[0] = mxCreateNumericArray(x0.getNdim(), x0.getSizes(), mxDOUBLE_CLASS, mxREAL);
    NdimData x(x0.getNdim(), x0.getSizes(), mxGetPr(plhs[0]));
    memcpy(x.getPtr(), x0.getPtr(), x0.getNumel() * sizeof(double));
    
    
    
    
    /////////////////////////
    //     computation     //
    /////////////////////////
    if ( (x.getSize(0) != AtY.getSize(0))  ||  (x.getSize(1) != AtY.getSize(1))  ||  (x.getNumSlices(*waveDims) != AtA.getSize(0))  ||  (x.getNumSlices(*waveDims) != AtA.getSize(1)) ) {
        error = true;
        RW_PRINT("ERROR: Incompatible sizes between S0, AtA and AtY\n");
    }
    if (!error)
        W.analysisInversion(x, AtA, AtY, thresholds, maximumIteration, nonNegative);
    
    
    if (waveDims)
        delete waveDims;
    if (filter)
        delete filter;

}










