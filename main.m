clc;
clear;
close all;


prm = getPrm(0);
        
folderName = 'res';
if exist(folderName,'dir')
    rmdir(folderName,'s');
end
mkdir(folderName);

index = 1;

for sigma = prm.noiseSigmas
    disp(sigma);
    for i = 1:prm.numOfImages
        disp(i);
        str = sprintf('test%2.2f\\N%2.2fIM%d.mat',sigma,sigma,i);
        load(str);
        im = Image(I);

        tic;
        im = im.orientedMeans;
        im = im.detectEdges(0,sigma);
        im = im.nonMaximalSupression(true);
        R = im.edgesImage;
        toc;

        figure,imshow(R>0);
        %figure,imshow(B);
        figure,imshow(I);
        
        cd res
        s.R = R;
        s.B = {B>0};
        s.I = I;
        
        save(sprintf('%d.mat',index), '-struct', 's');
        imwrite(R>0,sprintf('%dE.png',index),'PNG');
        imwrite(I,sprintf('%dI.png',index),'PNG');
        cd ..;
        index = index+1;
    end
end