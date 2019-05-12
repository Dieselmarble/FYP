%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	CDMA based using multiple PN sequences embeded into whole object
%           Watermark Recovery

clear all;

% save start time
start_time=cputime;

% read in the watermarked object
file_name='cdma_watermarked.bmp';
watermarked_image=double(imread(file_name));

% determine size of watermarked image
Mw=size(watermarked_image,1);	        %Height
Nw=size(watermarked_image,2);	        %Width

% read in original watermark
file_name='_copyright_small.bmp';
orig_watermark=double(imread(file_name));

% determine size of original watermark
Mo=size(orig_watermark,1);	%Height
No=size(orig_watermark,2);	%Width

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

% initalize message to all ones
message_vector=ones(1,Mo*No);
for kk=1:length(message_vector)

    % generate {-1,0,1} PN sequence    
    pn_sequence=round(2*(rand(Mw,Nw)-0.5));
    
    % calculate correlation
    if (watermarked_image == pn_sequence)
        correlation=1;
    else
        correlation(kk)=corr2(watermarked_image,pn_sequence);
    end
end
    
% use the average correlation value as threshold
threshold=mean(correlation);

% if correlation exceeds threshold, set message bit low
for kk=1:length(message_vector)
    if correlation(kk) > threshold
        message_vector(kk)=0;
    end
end

% reshape the message vector and display recovered watermark.
figure(2)
message=reshape(message_vector,Mo,No);
imshow(message,[])
title('Recovered Watermark')

% display processing time
elapsed_time=cputime-start_time,
