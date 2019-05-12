%Name:		Chris Shoemaker
%Course:	EER-280 - Digital Watermarking
%Project: 	Least Significant Bit Substitution
%           Watermark Embeding

clear all;

% save start time
start_time=cputime;

% read in the cover object
file_name='_lena_std_bw.bmp';
cover_object=imread(file_name);

% read in the message image
file_name='_copyright.bmp';
message=imread(file_name);

% convert to double for normalization, then back again
message=double(message);
message=round(message./256);
message=uint8(message);

% determine size of cover object
Mc=size(cover_object,1);	%Height
Nc=size(cover_object,2);	%Width

% determine size of message object
Mm=size(message,1);	        %Height
Nm=size(message,2);	        %Width

% title the message object out to cover object size to generate watermark
for ii = 1:Mc
    for jj = 1:Nc
        watermark(ii,jj)=message(mod(ii,Mm)+1,mod(jj,Nm)+1);
    end
end

% now we set the lsb of cover_object(ii,jj) to the value of watermark(ii,jj)
watermarked_image=cover_object;
for ii = 1:Mc
    for jj = 1:Nc
        watermarked_image(ii,jj)=bitset(watermarked_image(ii,jj),1,watermark(ii,jj));
    end
end

% write the watermarked image out to a file
imwrite(watermarked_image,'lsb_watermarked.bmp','bmp');

% display processing time
elapsed_time=cputime-start_time,

% calculate the PSNR
psnr=psnr(cover_object,watermarked_image,Mc,Nc),

% display watermarked image
figure(1)
imshow(watermarked_image,[])
title('Watermarked Image')