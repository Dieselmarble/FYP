close all;
% -------- add attacks --------------
% ------- gaussian noise attack ----------
A_ = imnoise(A_embed,'gaussian',0,0.03);
% ------- cropping attack ---------
% A_ = imcrop(A,[0 0 240 240]);
% A_ = padarray(A_,[240 240],0,'post');
% ------- gaussian filter ----------
% A_ = fspecial('gaussian',3,1);
% attackf=filter2(h,f);
% ------- salt & pepper noise -------
% A_ = imnoise(A,'salt & pepper',0.0002);
% reshape the image into patches
S_CI_=DvdBptSubBp(A_,k,l);
%--------------------------------------------------------------------------
icasig_=W_matrix*S_CI_;
index_=max_cov(icasig_);
Y_=[icasig_(index_,:); y2];

water = 2*Y_(2,:)-Y_(1,:);
water = normalize(water*2);
WWW = round(water,0);
figure;
imshow(reshape(WWW,64,64));

%--------------------------------------------------------------------------
[Watermark] = fastica(Y_);
Watermark = abs(Watermark);
% Wartermark = normalize(Wartermark);
WWW=ones(k,l);             
if mean(Watermark(2,:))<mean(Watermark(1,:))
    mmm=Watermark(2,:);
else mmm=Watermark(1,:); %choose the minimum row which contains the watermark
end
%--------------------------------------------------------------------------
for i=1:l
    WWW(:,i)=(mmm((i-1)*k+1:i*k))'; %0:k ; k+1:2k ; 2k+1 : 3k
end
%--------------------------------------------------------------------------
% remove negative entries
for i=1:l
    for j=1:k
        if WWW(j,i)<1    
            WWW(j,i)=1;
        else
            WWW(j,i)=0;
        end
    end
end
WWW_=255*(1-WWW);
WWW_ = 1- WWW;
figure;
imshow(WWW);
ratio = 10*log10(psnr(mark,WWW_,k,l));
ratio_ = 10*log10(psnr(mark,WWW,k,l));