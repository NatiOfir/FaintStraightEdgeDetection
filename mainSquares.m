clc;
clear;
close all;

param = 1;
prm = getPrm(param);
        
folderName = 'resSquares';
if exist(folderName,'dir')
    rmdir(folderName,'s');
end
mkdir(folderName);

index = 1;

sizes = [20 40 60];
snrs = 0.8:0.2:2;
noiseSigma = 0.1;

for i = sizes
    for j = snrs
        if j == 1 | j == 2
            src = sprintf('Sqr_4test\\src\\CC_4_SQR_%d_SNR_%dSrc+noise.png',i,j);
            gt = sprintf('Sqr_4test\\gt\\CC_4_SQR_%d_SNR_%dgt.png',i,j);
        else
            src = sprintf('Sqr_4test\\src\\CC_4_SQR_%d_SNR_%1.1fSrc+noise.png',i,j);
            gt = sprintf('Sqr_4test\\gt\\CC_4_SQR_%d_SNR_%1.1fgt.png',i,j);
        end
        I = im2double(imread(src));
        im = Image(I);
        B = imread(gt);
        tic;
        im = im.orientedMeans;
        signalSigma = noiseSigma*j;
        im = im.detectEdges(0,noiseSigma,signalSigma);
        im = im.nonMaximalSupression(true);
        %im = im.NMS(false);
        %im = im.createEdgesImage;
        R = im.edgesImage;
        toc;
        Topt = 0.01;
        figure,imshow(R>=Topt);
        figure,imshow(B);
        figure,imshow(I);
        
        cd resSquares
        s.R = R;
        s.B = {B};
        s.I = I;
        
        save(sprintf('%d.mat',index), '-struct', 's');
        imwrite(R>=Topt,sprintf('%dE.png',index),'PNG');
        imwrite(I,sprintf('%dI.png',index),'PNG');
        cd ..;
        index = index+1;
    end
end