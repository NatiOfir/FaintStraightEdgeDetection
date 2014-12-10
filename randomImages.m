function randomImages()
    clc;
    close all;
    
    prm = getPrm(1);
    
    for stdNoise = prm.noiseSigmas

        folderName = sprintf('%s%2.2f',prm.folderName,stdNoise);

        if exist(folderName,'dir')
            rmdir(folderName,'s');
        end
        mkdir(folderName);

        for i=1:prm.numOfImages
            I = zeros(prm.imSize)+prm.imMean;
            N = randn(prm.imSize)*stdNoise;
            I = I+N;
            B = zeros(prm.imSize);

            for j = 1:prm.numOfLines
                coords = randi([prm.imBorder prm.imSize-prm.imBorder],1,4);
                w = randi([1 prm.maxWidth],1);
                contrast = randn(1)*prm.stdContrast+prm.meanContrast;
                sign = randi([1 2],1)-2;
                if (sign == 0)
                    sign = 1;
                end
                contrast = contrast*sign;

                l = Line(coords(1),coords(2),coords(3),coords(4),w);
                l = l.calcPixels;
                L = l.getLineImage(prm.imSize,prm.imSize);
                b = l.getEdgeImage(prm.imSize,prm.imSize);
                B = B+b;
                L = L*contrast;
                I = I+L;    
            end

            s.I = I;
            s.B = B;
            s.noiseLevel = stdNoise;
            save(sprintf('%s/N%2.2fIM%d.mat',folderName,stdNoise,i), '-struct', 's');
            imwrite(I,sprintf('%s/N%2.2fIM%d.png',folderName,stdNoise,i),'PNG');

        end
    end
end

