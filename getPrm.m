function prm = getPrm(ind)

    prm.imSize = 200;
    prm.imBorder = 10;
    prm.numOfLines = 10;
    prm.numOfImages = 10;
    prm.maxWidth = 4;
    prm.meanContrast = 0.2;
    prm.stdContrast = 0.1;
    prm.imMean = 0.5;
    prm.folderName = 'test';
    prm.noiseSigmas = [0.01 0.1];
    
    switch ind
        case 0 %sines
            prm.mapFactor = 0.5;
            prm.p = 1.2;
            prm.s = 0.1;
            prm.sigmaNoise = 0.1;
            prm.useNormalSignal = false;
            prm.normalGap = 0.02;
            prm.normalA = 0.118;
            prm.normalB = 0.122;
            prm.useGaussSignal = false;
            prm.useBack = false;
            prm.doVarTest = false;
            prm.binaryResult = false;
            prm.useOldT = true;
            prm.onGaussSigma = 0.1;
            prm.pBack = 0.8;
            prm.sBack = 0.06;
            prm.alpha = 0.01;
            prm.filter = [1 0 -1]*1/2;
            prm.varTest = 0.01;
            prm.minContrast = 8; % power of 2, <=1 means no min Contrast
            prm.edgeSupressT = 100;%0.5; % the precentage of thersholding edges in the non maximal supression process. 
        case 1 %sines
            prm.mapFactor = 0.5;
            prm.p = 1.2;
            prm.s = 0.1;
            prm.useGaussSignal = false;
            prm.useBack = false;
            prm.doVarTest = false;
            prm.binaryResult = false;
            prm.onGaussSigma = 0.1;
            prm.pBack = 0.8;
            prm.sBack = 0.06;
            prm.alpha = 0.01;
            prm.filter = [1 1 0 -1 -1 ]*1/4;
            prm.varTest = 0.01;
            prm.minContrast = 1; % power of 2, <=1 means no min Contrast
            prm.edgeSupressT = 0.8; % the precentage of thersholding edges in the non maximal supression process. 
        case 2 % squares or sines
            prm.mapFactor = 0.5;
            prm.p = 1.2;
            prm.s = 0.1;
            prm.useGaussSignal = false;
            prm.useBack = false;
            prm.doVarTest = false;
            prm.binaryResult = false;
            prm.onGaussSigma = 0.1;
            prm.pBack = 0.8;
            prm.sBack = 0.06;
            prm.alpha = 0.01;
            %prm.filter = [0.5 0 -0.5];
            prm.filter = [1 1 0 -1 -1 ]*1/4;
            prm.varTest = 0.01;
            prm.minContrast = 8; % power of 2, <=1 means no min Contrast
            prm.edgeSupressT = 0.8; % the precentage of thersholding edges in the non maximal supression process. 
        otherwise
    end
end

