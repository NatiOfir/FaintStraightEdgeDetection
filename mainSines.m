clc;
clear;
close all;
        
folderName = 'resSines';
if exist(folderName,'dir')
    rmdir(folderName,'s');
end
mkdir(folderName);

index = 1;

sizes = [3 5 7];
snrs = 0.8:0.2:2;
noiseSigma = 0.1;

for i = sizes
    for j = snrs
        if j == 1 || j == 2
            src = sprintf('sin_4test\\src\\Sin_4_B_3_W_%d_SNR_%dSrc+noise.png',i,j);
            gt = sprintf('sin_4test\\gt\\Sin_4_B_3_W_%d_SNR_%dgt.png',i,j);
        else
            src = sprintf('sin_4test\\src\\Sin_4_B_3_W_%d_SNR_%1.1fSrc+noise.png',i,j);
            gt = sprintf('sin_4test\\gt\\Sin_4_B_3_W_%d_SNR_%1.1fgt.png',i,j);
        end
        I = im2double(imread(src));
        im = Image(I);
        B = imread(gt);
        tic;
        im = im.orientedMeans;
        signalSigma = noiseSigma*1.1;
        im = im.detectEdges(0,noiseSigma,signalSigma);
        im = im.nonMaximalSupression(true);
        %im = im.createEdgesImage;
        R = im.edgesImage;
        toc;
        Topt = 0.02;
        figure,imshow(R>=Topt);
        figure,imshow(B);
        figure,imshow(I);
        
        cd resSines
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