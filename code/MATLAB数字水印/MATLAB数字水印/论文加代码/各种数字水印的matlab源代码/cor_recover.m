%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Threshold-Based Correlation using block processing in the spatial domain
%           Watermark Recovery

clear all;

% save start time
start_time=cputime;

blocksize=16;       % set the size of the block in cover to be used for each bit in watermark

% read in the watermarked object
file_name='cor_watermarked.bmp';
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

% generate the watermark equal to the size of one block
pn_sequence=round(2*(rand(blocksize,blocksize)-0.5));

% pad message out to maximum message size with ones
message_vector=ones(No*Mo,1);

% process the image in blocks
% for each block determine it's correlation with base pn sequence
x=1;
y=1;
for (kk = 1:length(message_vector))

    % sets correlation to 1 when patterns are identical to avoid /0 errors
    % otherwise calcluate correlation 
    if (watermarked_image(y:y+blocksize-1,x:x+blocksize-1) == pn_sequence)
        correlation(kk)=1;
    else
        correlation(kk)=corr2(watermarked_image(y:y+blocksize-1,x:x+blocksize-1),pn_sequence);
    end
    
    % move on to next block. At and of row move to next row
    if (x+blocksize) >= Nw
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
end

% if correlation exceeds average correlation
for kk = 1:length(correlation)
    
    if (correlation(kk) > mean(correlation(1:Mo*No)))
        message_vector(kk)=0;
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