%%%%%%%%%%%%%%% Barbara image %%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;
imgt = double(imread('texture4.tif'));
imgc = double(imread('boy.tif'));
vecOrig1 = imgt(:)';
vecOrig2 = imgc(:)';



% normalise original image to unity variance
vecOrig1 = normalize(vecOrig1);
vecOrig2 = normalize(vecOrig2);

imgt = reshape(vecOrig1,[256,256]);
imgc = reshape(vecOrig2,[256,256]);

% Dictionary stuff (here Curvelets + UDWT).
dict1='CURVWRAP';pars11=2;pars12=0;pars13=0;
dict2='LDCT2';pars21=256;pars22=16/256;pars23=0; % Remove Low-frequencies < 16/256 from textured part.
dicts=MakeList(dict1,dict2);
pars1=MakeList(pars11,pars21);
pars2=MakeList(pars12,pars22);
pars3=MakeList(pars13,pars23);


% Call the MCA.
itermax 	= 100;
tvregparam 	= 0.1;
tvcomponent	= 2;
expdecrease	= 1;
lambdastop	= 1;
sigma		= 0;%1e-6;
display		= 0;

vecOrig1 = imgt(:)';
vecOrig2 = imgc(:)';
Am = rand(2,2);
x = Am*[vecOrig1;vecOrig2];

noise = randn(size(x));

PSNR_list=[3:5:20];
R_all = [];

for iter = 1:length(PSNR_list)
    PSNR = PSNR_list(iter);
    sigma = std(x(:))*10^(-PSNR/20);
    img = x + sigma*noise; 
    [parts,options]=mMCA2_Bcr(img,dicts,pars1,pars2,pars3,itermax,tvregparam,tvcomponent,expdecrease,lambdastop,[],sigma,display,2);
    options.inputdata = 'Input image: Boy + Texture 256 x 256';
    options
    [ST,I] = dbstack;
    name=eval(['which(''' ST(1).name ''')']);
    eval(sprintf('save %s options -V6',[name(1:end-2) 'metadata']));
    R1 = corr2(imgc,squeeze(parts(:,:,1)));
    R_all = [R_all;R1];
end

% % Display results.
% figure;
% set(gcf,'Name','MCA Texture + Boy','NumberTitle','off');
% subplot(321);
% imagesc(img);axis image;rmaxis;
% title('Original Texture + Boy');
% 
% subplot(322);
% imagesc(squeeze(sum(parts,3)));axis image;rmaxis;
% title(sprintf('Recovered MCA PSNR=%.3g dB',psnr(img,squeeze(sum(parts,3)))));
% 
% subplot(323);
% imagesc(imgc);axis image;rmaxis;
% title('Original Cartoon part');
% 
% subplot(324);
% imagesc(squeeze(parts(:,:,1)));axis image;rmaxis;
% title(sprintf('MCA Cartoon'));
% 
% subplot(325);
% imagesc(imgt);axis image;rmaxis;
% title('Original Texture');
% 
% subplot(326);
% imagesc(squeeze(parts(:,:,2)));axis image;rmaxis;
% title(sprintf('MCA Texture'));
% colormap('gray');
% 
% R1 = corr2(imgc,squeeze(parts(:,:,1)));
% 
