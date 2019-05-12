%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Comparison-Based Correlation using block processing in the spatial domain
%           Uses two PN sequences; one for a "0" and another for a "1"
%           Watermark Recovery

clear all;

% save start time
start_time=cputime;

blocksize=16;      % set the size of the block in cover to be used for each bit in watermark

% read in the watermarked object
file_name='cor_watermarked_mod.bmp';
watermarked_image=double(imread(file_name));

% determine size of watermarked image
Mw=size(watermarked_image,1);	%Height
Nw=size(watermarked_image,2);	%Width

% determine maximum possible message size in object
max_message=Mw*Nw/(blocksize^2);

% read in original watermark
file_name='_copyright.bmp';
orig_watermark=double(imread(file_name));

% determine size of original watermark
Mo=size(orig_watermark,1);	%Height
No=size(orig_watermark,2);	%Width

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

% generate PN sequences to designate "1" and "0"
watermark_one=round(2*(rand(blocksize,blocksize)-0.5));
watermark_zero=round(2*(rand(blocksize,blocksize)-0.5));

% find two highly un-correlated PN sequences
while (corr2(watermark_one,watermark_zero) > -0.1)
    watermark_one=round(2*(rand(blocksize,blocksize)-0.5));
    watermark_zero=round(2*(rand(blocksize,blocksize)-0.5));
end

% pad message out to maximum message size with ones
message_vector=ones(max_message,1);

% process the image in blocks
% for each block determine it's correlation with base pn sequence
x=1;
y=1;
for (kk = 1:length(message_vector))

    % calculate correlations for both PN sequences
    correlation_one(kk)=corr2(watermarked_image(y:y+blocksize-1,x:x+blocksize-1),watermark_one);
    correlation_zero(kk)=corr2(watermarked_image(y:y+blocksize-1,x:x+blocksize-1),watermark_zero);
    
    % choose which ever correlation is higher
    if correlation_one(kk) > correlation_zero(kk)
        message_vector(kk)=1;
    else
        message_vector(kk)=0;
    end
    
    % move on to next block. At and of row move to next row
    if (x+blocksize) >= Nw
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
end

% reshape the message
message=reshape(message_vector(1:Mo*No),Mo,No);

% display processing time
elapsed_time=cputime-start_time,

% display the recovered message
figure(2)
imshow(message,[])
title('Recovered Message')