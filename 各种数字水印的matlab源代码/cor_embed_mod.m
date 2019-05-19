%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Comparison-Based Correlation using block processing in the spatial domain
%           Uses two PN sequences; one for a "0" and another for a "1"
%           Watermark embedding

clear all;

% save start time
start_time=cputime;

k=5;            % set the gain factor for embeding
blocksize=16;   % set the size of the block in cover to be used for each bit in watermark

% read in the cover object
file_name='_lena_std_bw.bmp';
cover_object=double(imread(file_name));

% determine size of watermarked image
Mc=size(cover_object,1);	%Height
Nc=size(cover_object,2);	%Width

% determine maximum message size based on cover object, and blocksize
max_message=Mc*Nc/(blocksize^2);

% prefilter the cover image with a matched filter to reduce detection errors
%prefilt=[-1,-1,-1;-1,10,-1;-1,-1,-1]./2;
%cover_object=filter2(prefilt,cover_object);

% read in the message image
file_name='_copyright.bmp';
message=double(imread(file_name));
Mm=size(message,1);	%Height
Nm=size(message,2);	%Width

% reshape the message to a vector
message=round(reshape(message,Mm*Nm,1)./256);

% check that the message isn't too large for cover
if (length(message) > max_message)
    error('Message too large to fit in Cover Object')
end

% pad the message out to the maximum message size with 0's
message_vector=ones(1,max_message);
message_vector(1:length(message))=message;

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

% generate PN sequences to designate "1" and "0"
pn_sequence_one=round(2*(rand(blocksize,blocksize)-0.5));
pn_sequence_zero=round(2*(rand(blocksize,blocksize)-0.5));

% find two highly un-correlated PN sequences
while (corr2(pn_sequence_one,pn_sequence_zero) > -0.1)
    pn_sequence_one=round(2*(rand(blocksize,blocksize)-0.5));
    pn_sequence_zero=round(2*(rand(blocksize,blocksize)-0.5));
end

% process the image in blocks
% first construct the global watermark mask
x=1;
y=1;
for (kk = 1:length(message_vector))    
    
    % if message bit contains zero, add PN sequence to that portion of mask
    if (message_vector(kk) == 0)
        watermark_mask(y:y+blocksize-1,x:x+blocksize-1) = pn_sequence_zero;
        
    % otherwise mask is filled with zeros
    else
        watermark_mask(y:y+blocksize-1,x:x+blocksize-1) = pn_sequence_one;
    end
    
    % move to next block of mask along x; If at end of row, move to next row
    if (x+blocksize) >= Nc
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
end

% add watermark mask to cover image using gain factor k 
watermarked_image_dbl=cover_object+k*watermark_mask;
watermarked_image_int=uint8(cover_object+k*watermark_mask);

% write the watermarked image out to a file
imwrite(watermarked_image_int,'cor_watermarked_mod.bmp','bmp');

% display processing time
elapsed_time=cputime-start_time,

% calculate the PSNR
psnr=psnr(cover_object,watermarked_image_dbl,Mc,Nc),

% display watermarked image
figure(1)
imshow(watermarked_image_int,[])
title('Watermarked Image')

