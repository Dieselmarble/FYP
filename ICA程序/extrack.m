%--------------------------------------------------------------------------
% add attacks
close all;
A_ = imnoise(A,'gaussian',0,0.01);

A_ = fspecial('gaussian',3,1);
attackf=filter2(h,f);
att='gaussian low pass';

S_CI_=DvdBptSubBp(A_,k,l);
%--------------------------------------------------------------------------
icasig_=W_matrix*S_CI_;
index_=max_cov(icasig_);
Y_=[icasig_(index_,:); y2];
%--------------------------------------------------------------------------
[Wartermark] = fastica(Y_);
%imshow(Wartermark);
Wartermark=abs(Wartermark);
WWW=ones(k,l);             
if mean(Wartermark(2,:))<mean(Wartermark(1,:))
    mmm=Wartermark(2,:);
else mmm=Wartermark(1,:); %choose the minimum row
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
figure;
imshow(WWW);
ratio = psnr(mark,WWW,k,l);