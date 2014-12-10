function [ imgCell,curvesCell ] = readAllData()
    % this function reads all BSD data into 2 cells
    % imgCell = cell containning all the 500 images
    % curvesCell = cell containning the sum of all the human bounderies
    dirData = dir('BSR\BSDS500\data\images\all');
    fileList = {dirData.name};
    fileList(1:2) = [];
    numOfImages = numel(fileList);
    imgCell = cell(1,numOfImages);
    curvesCell = cell(1,numOfImages);

    for i=1:numOfImages
        if(mod(i,numOfImages/50) == 0)
            disp(sprintf('Image = %d',i));
        end
        cd BSR\BSDS500\data\images\all;
        imgCell{i} = rgb2gray(im2double(imread(fileList{i})));
        cd ..\..;
        cd groundTruth\all;
        spl = regexp(fileList{i}, '\.', 'split');
        load(strcat(spl{1},'.mat'));
        cd ..\..\..\..\..;
        curvesCell{i} = cell(1,0);
        for j=1:numel(groundTruth)
            curvesCell{i}{j} = groundTruth{j}.Boundaries;
        end
    end
    
    s.imgCell = imgCell;
    s.curvesCell = curvesCell;
    
    save('data.mat', '-struct', 's');
end

