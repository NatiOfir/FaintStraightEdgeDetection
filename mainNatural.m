clc;
clear;
close all;

load('data.mat');
N = 10;
sigma = 0.05;

folderName = 'resNatural';
if exist(folderName,'dir')
    rmdir(folderName,'s');
end
mkdir(folderName);

index = 1;


for i=1:1
    I = Image.normalize(imgCell{i});
    B = curvesCell{i};
    im = Image(imgCell{i});

    tic;
    im = im.orientedMeans;
    im = im.detectEdges(0,sigma,sigma);
    size(im.edges)
    im = im.nonMaximalSupression(true);
    %im = im.createEdgesImage;
    R = im.edgesImage;
    toc;

    sumB = zeros(size(B{1}));
    for j = 1:numel(B)
        sumB = sumB+B{j};
    end
    
    figure,imshow(R>0);
    figure,imshow(sumB);
    figure,imshow(I);
    
    cd resNatural
    s.R = R;
    s.B = B;
    s.I = I;
    

    save(sprintf('%d.mat',index), '-struct', 's');
    imwrite(R>0,sprintf('%dE.png',index),'PNG');
    imwrite(I,sprintf('%dI.png',index),'PNG');
    imwrite(sumB,sprintf('%dB.png',index),'PNG');

    cd ..;
    index = index+1;
    
end


