% initialization of the nGMCAlab toolbox

if ~exist('nGMCA_alg', 'file')
    
    curFolder = [pwd() '/'];
    getd = @(rel_p) path([curFolder rel_p], path);
    
    %main folder
    getd('');
    
    %useful tools which are commonly used
    getd('Tools')
    
    %proximal algorithms and methods
    getd('ProximalMethods')
    
    %NMF algorithms
    getd('NMF')
    getd('NMF/nGMCA/bricks')
    getd('NMF/nGMCA')
    getd('NMF/nGMCA/test_versions') %tests with hard thresholding and naive implementation
    getd('NMF/NMFtools')
    
    %used to perform benchmarks
    getd('Benchmarks')
    getd('Benchmarks/BenchmarkTools')
    
    %external libraries
    getd('Libraries/munkres')%source assignment for evaluation
    
    %wavelets
    getd('Libraries/redWaveToolbox/bin')

    
end

if ~exist('MakeONFilter', 'file')
    fprintf(1, 'This toolbox requires the WaveLab toolbox to work correctly\n');
    fprintf(1, 'when using wavelet transforms. Please consider downloading it\n');
    fprintf(1, 'and adding it to your path.\n');
    fprintf(1, 'More specifically, it uses its implementation of the orthogonal\n');
    fprintf(1, 'wavelets and the MakeONFilter function to set the wavelet filters.\n');
end

% unify random number generator seed
if ~exist('rng', 'var')
    randGenSeed = @(i) rand('twister', i);
else
    randGenSeed = @(i) rng(i, 'twister');
end

clear getd curFolder;
