function randomLines()
    clc;
    close all;    
    folderName = 'Line_4test';

    if exist(folderName,'dir')
        rmdir(folderName,'s');
    end
    mkdir(folderName);
    size = 129;
    numOfLines = 1;
    coords = randi([15 129-15],numOfLines,4);
    prm = getPrm(0);
    
    for w = [4 6 8]; 
        for i=0.8:0.2:2
            I = zeros(size)+0.5;
            N = randn(size)*prm.sigmaNoise;
            I = I+N;
            B = zeros(size);
            
            for j = 1:numOfLines
                sign = randi(2,1)*2-3;
                contrast = (rand(1)*(prm.normalB - prm.normalA)+prm.normalA)*sign;
                l = Line(coords(j,1),coords(j,2),coords(j,3),coords(j,4),w);
                l = l.calcPixels;
                L = l.getLineImage(129,129);
                b = l.getEdgeImage(129,129);
                B = max(B,b);
                L = L*contrast;
                I = I+L;    
            end

            figure, imshow(B);
            figure,imshow(I);
            s.I = I;
            s.B = B;
            save(sprintf('%s/N%1.1fW%df.mat',folderName,i,w), '-struct', 's');
            imwrite(I,sprintf('%s/N%1.1fW%d.png',folderName,i,w),'PNG');
            imwrite(B,sprintf('%s/N%1.1fW%db.png',folderName,i,w),'PNG');
        end
    end
end

