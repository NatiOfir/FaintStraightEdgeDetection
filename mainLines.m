clc;
clear;
close all;

prm = getPrm(0);
        
folderName = 'resLines';
if exist(folderName,'dir')
    rmdir(folderName,'s');
end
mkdir(folderName);

index = 1;
noiseSigma = prm.sigmaNoise;

for w = [4 6 8]
    for i=0.8:0.2:2
        load(sprintf('%s/N%1.1fW%df.mat','Line_4test',i,w));

        I = im2double(I);
        im = Image(I);
        tic;
        im = im.orientedMeans;
        signalSigma = noiseSigma*1.1;
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
        
        cd resLines
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