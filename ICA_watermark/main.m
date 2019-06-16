%-------------------????------------------------------------------------
clear all; close all; clc;
timer = 0;
while timer < 1
c=0.3;
a=imread('lina.jpg'); %original image
b=imread('changsha.bmp')*255; %watermark
[m1,n1]=size(a);
[m2,n2]=size(b);hnjmk
e0=(sum(sum(a.^2)))/(m1*n1);
e0=c*e0;%intensity of watermark
[t,tkey]=sdwt(double(a),'db2',m2,n2,e0);%??????
[t,tkey]=embed(t,tkey,b);%????
aw=sidwt(t,'db2');%??
figure(1);
subplot(1,2,1);imshow(uint8(a));title('original image');
subplot(1,2,2);imshow(uint8(aw));title('embedded image');
imwrite(uint8(aw),'watermark.jpg');
% csvwrite('key.txt',reshape(tkey,m2,n2));
v1=m1*m1*255*255;
v2=sum(sum((double(a)-aw).^2));
snr=10*log10(v1/v2);% signal to noise ratio
disp('SNR');
disp(snr);
%---?????--------------------------------------------------------------
% clear;
% pause; 
f=imread('watermark.jpg'); 
%??????f???,????????
m=max(max(f));
f=double(f)./double(m);
%---??-------------------------------------------------------------------
attack=1;
switch attack
    case 0
        attackf=f;
        att='???';
    case 1    
%%1. JPEG ??
 imwrite(f,'attackf.jpg','jpg','quality',30);
 attackf=imread('attackf.jpg');
 attackf=double(attackf)/255;
 att='JPEG compress';
    case 2
% %2. ??????
h=fspecial('gaussian',3,1);
attackf=filter2(h,f);
att='gaussian low pass';
    case 3
%%3. ??????
attackf=histeq(f);
att='histogram average';
    case 4
%%4. ????
attackf=imadjust(f,[],[0.4,1]);
att='image lightening';
    case 5
%%5. ????
attackf=imadjust(f,[],[0,0.85]);
att='iamge dimming';
    case 6
%%6. ?????
attackf=imadjust(f,[0.3,0.6],[]);
att='increase contrast';
    case 7
%%7. ?????
attackf=imadjust(f,[],[0.2,0.8]);
att='decrease contrast';
    case 8
%%8. ??????
attackf=imnoise(f,'gaussian',0,0.01);
att='add gaussian noise';
    case 9
%%9. ????
attackf=imnoise(f,'salt & pepper',0.06);
att='add pepper noise';
    case 10
%%10. ???????
attackf=imnoise(f,'speckle',0.08);
att='add multiplicative noise';
end
%---?????--------------------------------------------------------------
f=attackf.*double(m);
figure(2);
imshow(uint8(f));%????????????
title(att);
imwrite(uint8(f),'watermark.jpg');
%---????----------------------------------------------------------------
% clear;
a=imread('watermark.jpg');
t=sdwt_ex(double(a),'db2',tkey);%???????
[w,map]=extract(t,tkey);%????
[r,c]=size(w);
figure(3);
for i=1:r
    subplot(ceil(r/3),3,i)
    imshow(255-100*abs(uint8(reshape(w(i,:),map(1),map(2)))));
    title(strcat('extract watermark',num2str(i)));
end;
%---???????----------------------------------------------------------
Rarr=[];
for i=1:r   %????????
    for j=i+1:r
       R(i,j)=abs(corr2(reshape(w(i,:),map(1),map(2)),reshape(w(j,:),map(1),map(2))));
       if i~=j
           Rarr=[Rarr;i,j,abs(R(i,j))];
       end;
    end;
end;

[cor,idx]=sort(Rarr(:,3),'descend');
max1=cor(1);max2=cor(2);
m1=Rarr(idx(1),1);n1=Rarr(idx(1),2);
m2=Rarr(idx(2),1);n2=Rarr(idx(2),2);
maxcor(1)=abs(corr2(reshape(w(m1,:),map(1),map(2)),reshape(w(m2,:),map(1),map(2))));
maxcor(2)=abs(corr2(reshape(w(m1,:),map(1),map(2)),reshape(w(n2,:),map(1),map(2))));
maxcor(3)=abs(corr2(reshape(w(n1,:),map(1),map(2)),reshape(w(m2,:),map(1),map(2))));
maxcor(4)=abs(corr2(reshape(w(n1,:),map(1),map(2)),reshape(w(n2,:),map(1),map(2))));
if mean(maxcor)>(max1*0.8)
    stdw=fuse_pca(reshape(w(m1,:),map(1),map(2)), reshape(w(n1,:),map(1),map(2)));
    figure(4);
    subplot(1,2,1)
    imshow(b);
    title('original watermark');
    subplot(1,2,2)
    wf=255-100*abs(uint8(stdw));
    imshow(wf);
    title('watermark after pca fusion');
else
    disp('fusion again');
end;
if abs(corr2(wf,b))>0.9
    break;
end;
    timer = timer + 1;
end;
disp('similarity between watermarks');
disp(corr2(wf,b));
%--------------------------------------------------------------------------