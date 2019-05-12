%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	CDMA based using multiple PN sequences embeded into whole object
%           Watermark Embeding

clear all;

% save start time
start_time=cputime;

k=2;                % set the gain factor for embeding

% read in the cover object
file_name='_lena_std_bw.bmp';
cover_object=double(imread(file_name));

% determine size of watermarked image
Mc=size(cover_object,1);	%Height
Nc=size(cover_object,2);	%Width

% read in the message image and reshape it into a vector
file_name='_copyright_small.bmp';
message=double(imread(file_name));
Mm=size(message,1);	                        %Height
Nm=size(message,2);	                        %Width
message_vector=round(reshape(message,Mm*Nm,1)./256);

% read in key for PN generator
file_name='_key.bmp';
key=double(imread(file_name))./256;

% reset MATLAB's PN generator to state "key"
rand('state',key);

watermarked_image=cover_object;

% when message contains a '0', add pn sequence with gain k to cover image
for kk=1:length(message_vector)
    pn_sequence=round(2*(rand(Mc,Nc)-0.5));
    
    if message(kk) == 0
        watermarked_image=watermarked_image+k*pn_sequence;
    end
end

% convert back to uint8
watermarked_image_uint8=uint8(watermarked_image);

% write watermarked Image to file
imwrite(watermarked_image_uint8,'cdma_watermarked.bmp','bmp');

% display processing time
elapsed_time=cputime-start_time,

% calculate the PSNR
psnr=psnr(cover_object,watermarked_image_uint8,Mc,Nc),

% display watermarked Image
figure(1)
imshow(watermarked_image_uint8,[])
title('Watermarked Image')
