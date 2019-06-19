close all; clear all; clc;

SNR_db = 30;  %-- SNR in dB
n = 128;
n0 = [];

iS = 256;
nTRs = 256;
nComps = 4;
nMixtures = 10;
pixles = 65536;
%sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im1 = double(imread('boat256.png'));
im2 = double(imread('paraty256_2.tif'));
im3 = double(imread('barbara256.png'));
im4 = double(imread('pakhawaj.tif'));

M0 = im3;
M1 = im2;
M0 = rescale(crop(M0,n)); 
M0 = rescale(sum(M0,3) );

M1 = rescale(crop(M1,n)); 
M1 = rescale(sum(M1,3) );

M = [M0,M1];

% options common to all the transforms
opt.Jmin = 3;   % minimum scale for wavelets
opt.n = n;      % size of the image

components = {};
%% load the wavelet dictionary
opt.threshold_factor = 1; % reduce influency if <1
clear cpt;
cpt.options = opt;
cpt.callback =  @callback_atrou;
cpt.name = 'wav';
cpt.tv_correction = 1; % add TV penalty
components{1} = cpt;
%% load the local DCT dictionary
opt.w = 32;      % size of patches
opt.q = opt.w/4; % controls overlap
opt.threshold_factor = 1; % reduce influency if <1
clear cpt;
opt.dct_type = 'orthogonal4';
opt.dct_type = 'redundant';
opt.remove_lowfreq = 1; % force to 0 the lowfrequency
cpt.options = opt;
cpt.callback =  @callback_localdct;
cpt.name = 'dct';
components{end} = cpt;

%% options for MCA
options.niter = 60;
options.Tmax = 2.5;
options.Tmin = 0;
options.n = n;
% options.InpaintingMask = U;

ML = perform_mca(M, components, options);