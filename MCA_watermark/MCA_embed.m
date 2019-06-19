%%%%%%%%%%%%%%% Barbara image %%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;
k=20;                           % set gain factor for embeding
blocksize=8;                    % set the dct blocksize

midband=[   0,0,0,1,1,1,1,0;    % defines the mid-band frequencies of an 8x8 dct
            0,0,1,1,1,1,0,0;
            0,1,1,1,1,0,0,0;
            1,1,1,1,0,0,0,0;
            1,1,1,0,0,0,0,0;
            1,1,0,0,0,0,0,0;
            1,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0 ];

% read the cover image
imgc = double(imread('lena.jpg'));
% read in the message image
file_name='changsha.bmp';
message=double(imread(file_name))*255;

Nw=size(message,1);	                %Height
Mw=size(message,2);	                %Width
Nc = size(imgc,1);	
Mc = size(imgc,2);
% normalise original image to unity variance

% determine maximum message size based on cover object, and blocksize
max_message=Mc*Nc/(blocksize^2);

% read in the message image
file_name='changsha.bmp';
message=double(imread(file_name))*255;
Mm=size(message,1);	                %Height
Nm=size(message,2);	                %Width

% reshape the message to a vector
message=round(reshape(message,Mm*Nm,1)./256);

% check that the message isn't too large for cover
if (length(message) > max_message)
    error('Message too large to fit in Cover Object')
end

% pad the message out to the maximum message size with ones's
message_vector=ones(1,max_message);
message_vector(1:length(message))=message;

% pad the message out to the maximum message size with ones's
message_vector=ones(1,max_message);
message_vector(1:length(message))=message;

%% Dictionary stuff (here Curvelets + UDWT).
dict1='CURVWRAP';
pars11=2;pars12=0;pars13=0;
dict2='LDCT2';
pars21=256;pars22=16/256;pars23=0; % Remove Low-frequencies < 16/256 from textured part.
dicts=MakeList(dict1,dict2);
pars1=MakeList(pars11,pars21);
pars2=MakeList(pars12,pars22);
pars3=MakeList(pars13,pars23);

%% Call the MCA.
itermax 	= 100;
tvregparam 	= 0.1;
tvcomponent	= 1;
expdecrease	= 1;
lambdastop	= 1;
sigma		= 0;
display		= 0;
% [parts,options]=MCA2_Bcr(imgc,dicts,pars1,pars2,pars3,itermax,tvregparam,tvcomponent,expdecrease,lambdastop,[],sigma,display);
% options.inputdata = 'Input image: Boy + Texture 256 x 256';
load('decompose.mat');
[ST,I] = dbstack;
% name=eval(['which(''' ST(1).name ''')']);
% eval(sprintf('save %s options -V6',[name(1:end-2) 'metadata']));
%%
% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;
key = sum(key);
% reset MATLAB's PN generator to state "key"
rand('state',key);

% generate PN sequence
pn_sequence_zero=round(2*(rand(1,sum(sum(midband)))-0.5));

% generate shell of watermarked image
cover_object = parts(:,:,1); % contour
watermarked_image = cover_object;
% process the image in blocks
x=1;
y=1;
for (kk = 1:length(message_vector))

    % transform block using DCT
    dct_block=dct2(cover_object(y:y+blocksize-1,x:x+blocksize-1));
    
    % if message bit contains zero then embed pn_sequence_zero into the mid-band
    % componants of the dct_block
    ll=1;
    if (message_vector(kk)==0)
        for ii=1:blocksize
            for jj=1:blocksize
                if (midband(jj,ii)==1) % if the resolution at this point is 1
                    dct_block(jj,ii)=dct_block(jj,ii)+k*pn_sequence_zero(ll);
                    ll=ll+1;
                end
            end
        end
    end
    
    % transform block back into spatial domain
    watermarked_image(y:y+blocksize-1,x:x+blocksize-1)=idct2(dct_block);    
    
    % move on to next block. At and of row move to next row
    if (x+blocksize) >= Nc
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
end
%%

% generate and display watermarked image
temp = watermarked_image + parts(:,:,2);
figure;
imshow(temp,[])
title('Watermarked Image')
% convert to uint8 and write the watermarked image out to a file
watermarked_image_int=uint8(temp);
imwrite(watermarked_image_int,'MCA_watermarked.bmp','bmp');
%
psnr_c=10*log10(psnr(imgc,watermarked_image_int,Nc,Mc));%%
