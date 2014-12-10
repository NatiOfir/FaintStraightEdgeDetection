function o=MarkBib2im(im,M)

o1=im2double(im);
o2=im2double(im);
o3=im2double(im);

o1(M==1)=1;
o2(M==1)=0;
o3(M==1)=0;
o(:,:,1)=o1;
o(:,:,2)=o2;
o(:,:,3)=o3;


