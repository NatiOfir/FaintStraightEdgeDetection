classdef Image
    % A class for a gray scale image
    % The class methods can detect all edges in an image and store them.
    % The line scanning is working in NlogN efficent algorithm.
    
    properties(GetAccess = public)
        % The grayscale image
        I
        
        % Edges x 7 matrix where every row is an edge data (x0,y0,x1,y1,score,length,isVertical)
        % sorted according to length firstly , and score as secondary.
        edges
        
        % A gray scale image which visualize the matrix edge, pixel value
        % is the maximum score it has in edges.
        edgesImage
        
        % A struct stores all the vertical oriented means of the image
        vertical
        
        % A struct stores all the horizontal oriented means of the image
        horizontal
        
        % The parameter index to use
        params
        
        % The variance of the assumed noise in the input image. i.i.d
        % addative gausssian noise is assumed.
        noiseSigma
        
        % The variance of the assumed signal in the input image.
        signalSigma
    end 
    properties(Constant)
        % width of a line, being used for edge drawing
        WIDTH = 1;
        EDGE_COLS = 7;
    end 
    
    methods (Access = public)
        % Constructor
        function obj = Image(I)
            obj.I = I;
            %obj.I = obj.normalize(I);
            obj.edges = cell(0);
        end
        
        % Calculates all oriented Means responses of the image. Stores them
        % in the structs vertical an horizantal.
        function obj = orientedMeans(obj)
            obj.vertical = obj.calcResponses(obj.I);
            obj.horizontal = obj.calcResponses(obj.I');
        end
        
        % Detect all edges in the images. Assumes that orientedMeans() was
        % called before.
        % params = the parametrs index of getPrm()
        % noiseSigma = the noise variance level in the image
        function obj = detectEdges(obj,params,noiseSigma,signalSigma)
            obj.params = params;
            obj.noiseSigma = noiseSigma;
            obj.signalSigma = signalSigma;
            obj.edges = zeros(0,obj.EDGE_COLS);
            obj = obj.detectEdgesFromResp(obj.vertical,true);
            obj = obj.detectEdgesFromResp(obj.horizontal,false);
            
            % sort edges table
            [Y,I1] = sort(obj.edges(:,obj.EDGE_COLS-2));
            lenCol = obj.edges(I1,obj.EDGE_COLS-1);
            [Y,I2] = sort(lenCol);
            
            obj.edges = obj.edges(I1,:);
            %obj.edges = obj.edges(I2,:);
            obj.edges = obj.edges(end:-1:1,:);
        end
        
        % Create the drawing of the edge scores. Assumes that detectEdges()
        % was called before.
        function obj = createEdgesImage(obj)
            prm = getPrm(obj.params);
            obj.edgesImage = zeros(size(obj.I));
            [m,n] = size(obj.edges);
            [M,N] = size(obj.I);
            for i=1:m
                l = Line(obj.edges(i,1),obj.edges(i,2) ,obj.edges(i,3),obj.edges(i,4),obj.WIDTH);
                l = l.calcPixels;
                obj.edgesImage = max(obj.edgesImage,l.getLineImage(M,N)*obj.edges(i,5));
            end
            
            if prm.binaryResult
                obj.edgesImage = obj.edgesImage > 0;
            else
                obj.edgesImage = obj.normalize(obj.edgesImage);
            end   
        end
        
        % calculates the onEdge probabilities of the input samples x.
        function prob = onEdgeProb(obj,x)
            prm = getPrm(obj.params);
            s = prm.s;
            p = prm.p;
            alpha = prm.alpha;
            z = 2*s/p*gamma(1/p)-alpha*(1+exp(-(alpha/s)^p));
            prob = -abs(x./s).^p-log(z);
            prob(abs(x)<=alpha) = -inf;
        end
        
        % calculates the onEdge probabilities of the input samples x.
        function prob = onEdgeGaussProb(obj,x,sigma)
            prob = log(normpdf(x,0,sigma));
        end
        
        function prob = onEdgeNormalProb(obj,x,a,b)
            a = a/2;
            b = b/2;
            inRange = (abs(x)<=b) & (abs(x) >=a);
            prob = log(inRange./(2*(b-a)));
        end
        
        % calculates the noise probabilities of input samples x of length len.
        function prob = noiseProb(obj,x,len)
            prm = getPrm(obj.params);
            w = sum(abs(prm.filter)>0);
            sigmaL = obj.noiseSigma^2/(w*len);
            prob = -0.5*log(2*pi*sigmaL)-x.^2/(2*sigmaL);
        end
        
        % calculate the offEdge probabilities of input samples x of length
        % len
        function prob = offEdgeProb(obj,x,len)
            prm = getPrm(obj.params);
            s = prm.sBack/len;
            p = prm.pBack;
            z = 2*s/p*gamma(1/p);
            prob = -abs(x./s).^p-log(z);
        end  
        
        % eliminate all the non maximal edges from the edge table
        function obj = nonMaximalSupression(obj,mex)
            prm = getPrm(obj.params);
            if mex
                [obj.edges obj.edgesImage] = NonMaximalSupression_mex(obj.edges,obj.I,prm.edgeSupressT);
            else
                [obj.edges obj.edgesImage] = NonMaximalSupression(obj.edges,obj.I,prm.edgeSupressT);
            end
            
            
            if prm.binaryResult
                obj.edgesImage = obj.edgesImage > 0;
            else
                obj.edgesImage = obj.normalize(obj.edgesImage);
            end
        end
        
        % Meirav NMS
        function obj = NMS(obj,mex)
            prm = getPrm(obj.params);
            EdgeList = obj.edges;
            isVertical = EdgeList(:,7);
            isVertical(isVertical ==0) = 3;
            EdgeList = [log2(EdgeList(:,6))-1 isVertical EdgeList(:,1) EdgeList(:,2) EdgeList(:,3) EdgeList(:,4) EdgeList(:,8)];
            EdgeListFinal = NMS(EdgeList,prm.filter,obj.noiseSigma,prm.edgeSupressT,obj.I);
            obj.edges = [EdgeListFinal(:,3) EdgeListFinal(:,4) EdgeListFinal(:,5) EdgeListFinal(:,6)];
        end
    end
    
    methods (Access = private)
        % detect the edges from the oriented means response. Assumes that
        % oriented means was called before.
        % Edge critera: onProb > offProb
        function obj = detectEdgesFromResp(obj,resp,isVertical)
            N = numel(obj.I);
            prm = getPrm(obj.params);
            for j = 1:length(resp)
                len = 2^(j-2);
                
                if(len<=2)
                    continue;
                end

                minContrast = resp{j}.minContrast;
                maxContrast = resp{j}.maxContrast;
                diffs = resp{j}.diffs;
                diffSquares = resp{j}.diffSquares;
                diffs = diffs/len;
                maxHalf = resp{j}.maxHalf;
                minHalf = resp{j}.minHalf;
                maxFourth = resp{j}.maxFourth;
                minFourth = resp{j}.minFourth;
                maxEigth = resp{j}.maxEigth;
                minEigth = resp{j}.minEigth;
                
                % compare probabilities
                if prm.useBack
                    highDiffs = obj.onEdgeProb(diffs) - obj.offEdgeProb(diffs,len);
                elseif prm.useGaussSignal
                    %highDiffs = obj.onEdgeGaussProb(diffs,obj.signalSigma) - obj.noiseProb(diffs,len);                    
                    highDiffs = obj.onEdgeGaussProb(diffs,prm.onGaussSigma) - obj.noiseProb(diffs,len); 
                elseif prm.useNormalSignal
                    a = obj.signalSigma-prm.normalGap;
                    b = obj.signalSigma+prm.normalGap;
                    %a = obj.signalSigma-2*obj.noiseSigma/len;
                    %b = obj.signalSigma+2*obj.noiseSigma/len;
                    highDiffs = obj.onEdgeNormalProb(diffs,a,b) - obj.noiseProb(diffs,len); 
                elseif prm.useOldT
                    w = sum(abs(prm.filter)>0);
                    highDiffs = abs(diffs) - obj.noiseSigma*sqrt(2*log(N)/(w*len));
                    %highDiffs = abs(diffs) - sqrt((obj.noiseSigma^2/(w*len-1))*log((N-1)*w*len));
                else
                    highDiffs = obj.onEdgeProb(diffs) - obj.noiseProb(diffs,len);
                end
                
                
                % min contrast test
                if prm.minContrast > 1
                    if len > prm.minContrast*2
                        minLen =  prm.minContrast;
                        minInRange = (minContrast/minLen >= diffs/2) & diffs>0;
                        maxInRange = (maxContrast/minLen <= diffs/2) & diffs<0;
                        minContrastInRange = minInRange | maxInRange;
                    else
                        varEst = diffSquares/(len-1) - diffs.^2*len/(len-1);
                        varInRange = abs(varEst./(obj.noiseSigma^2/2)-1) <= sqrt(2)*prm.varTest^(-0.5)/sqrt(len-1);
                        minContrastInRange = varInRange;
                    end
                else
                    if len >= 64
                        curLen = len/8;
                        minRange = (minEigth/curLen >= diffs/2) & diffs>0;
                        maxRange = (maxEigth/curLen <= diffs/2) & diffs<0;
                        minContrastInRange = minRange | maxRange;
                    elseif len >= 32
                        curLen = len/4;
                        minRange = (minFourth/curLen >= diffs/2) & diffs>0;
                        maxRange = (maxFourth/curLen <= diffs/2) & diffs<0;
                        minContrastInRange = minRange | maxRange;
                    elseif len >=4
                        varEst = diffSquares/(len-1) - diffs.^2*len/(len-1);
                        varInRange = abs(varEst./(obj.noiseSigma^2/2)-1) <= sqrt(2)*prm.varTest^(-0.5)/sqrt(len-1);
                        minContrastInRange = varInRange;
                    else
                        minContrastInRange = zeros(size(diffs));
                    end
                end
                
                if prm.useOldT
                    isEdge = (highDiffs >  0) & minContrastInRange; 
                else
                    isEdge = (highDiffs >  log(N^prm.mapFactor-1)) & minContrastInRange; 
                end
                
                [x,y,z] = ind2sub(size(highDiffs),find(isEdge));
                
                % read response values
                values = highDiffs(isEdge);
                %values = exp(values);
                
                edgesInd = obj.getXYIndices(isVertical,j,x,y,z,values);
                obj.edges = [obj.edges; edgesInd];
            end
        end
    end
    
    methods(Static)
        % convert a index an value in the oriented means to x,y points and
        % edge tuple. 
        function ind = getXYIndices(isVertical,cellIndex,curX,curY,angle,values)
            length = floor(2^(cellIndex-2));
            lenVec = zeros(size(values))+length;

            verticalVec = zeros(size(values))+isVertical;
            
            if length == 0
                x0 = curX;
                x1 = curX;
                y0 = curY;
                y1 = curY;
            else
                if length>2
                    vStep = 2;
                else
                    vStep = 1;
                end
                y0 = curY;
                x0 = 1+(curX-1).*length/vStep;
                x1 = x0+length;
                center = length+1;
                rel = angle-center;
                y1 = y0+rel;
            end
            
            if isVertical
               ind = [x0 y0 x1 y1 values lenVec verticalVec];
            else
                ind = [y0 x0 y1 x1 values lenVec verticalVec];
            end
        end
        
        % calculates all the vertical responses of an image in NlogN
        % approximation.
        function resp = calcResponses(img)
            prm = getPrm(0);
            [m,n] = size(img);
            iter = floor(log2(m-1))+1;
            %iter = 7;
            resp = cell(1,iter+1);
            s.sums = 0.5*img;
            s.squares = 0.5*(img.^2);
            s.diffs = imfilter(img,prm.filter,'symmetric');
            s.diffSquares = 0.5*s.diffs.^2;
            s.diffs = 0.5*s.diffs;
            s.minContrast = s.diffs;
            s.maxContrast = s.diffs;
            
            resp{1} = s;
            
            for i=1:iter
                s = [];
                length = 2^(i-1);
                % number of angles of that lenght
                angles = length*2+1;
                % center angle
                center = length+1;
                %last layer responses
                prev = resp{i};
                [m,n,z] = size(prev.sums);
                
                % set the vertical step between samples. if length>=4, we
                % jump 2 lines at a time.
                if length>2
                    vStep = 2;
                else
                    vStep = 1;
                end
                
                rows = floor((m-1)/vStep);
                s.sums = zeros(rows,n,angles);
                s.squares = zeros(rows,n,angles);
                s.diffs = zeros(rows,n,angles);
                s.diffSquares = zeros(rows,n,angles);
                s.minContrast = zeros(rows,n,angles);
                s.maxContrast = zeros(rows,n,angles);
                s.minHalf = zeros(rows,n,angles);
                s.minFourth = zeros(rows,n,angles);
                s.maxHalf = zeros(rows,n,angles);
                s.maxFourth = zeros(rows,n,angles);
                s.minEigth = zeros(rows,n,angles);
                s.maxEigth = zeros(rows,n,angles);
                
                % indices for the angles. j1 indicates the angle of the
                % upper half. j2 indicates the angle of the lower half.
                j1=1;
                j2=1;
                for j = 1:angles
                    % the angle, relative to the center angle.
                    rel = j-center;
                    % padding for the boundaries.
                    pad = zeros(rows,abs(rel),1);
                    %the relative angle of the upper half.
                    rel2 = ceil(abs(rel)/2);
                    
                    % angle to the small, padding from the small boundary.
                    if rel < 0
                        rel = -rel;
                        s.sums(:,:,j) = [pad prev.sums(1:vStep:end-vStep,1+rel:end,j1)+prev.sums(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2)];
                        s.squares(:,:,j) = [pad prev.squares(1:vStep:end-vStep,1+rel:end,j1)+prev.squares(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2)];
                        s.diffs(:,:,j) = [pad prev.diffs(1:vStep:end-vStep,1+rel:end,j1)+prev.diffs(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2)];
                        s.diffSquares(:,:,j) = [pad prev.diffSquares(1:vStep:end-vStep,1+rel:end,j1)+prev.diffSquares(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2)];
                        
                        if length <= prm.minContrast
                            s.minContrast(:,:,j) = [pad prev.diffs(1:vStep:end-vStep,1+rel:end,j1)+prev.diffs(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2)];
                            s.maxContrast(:,:,j) = s.minContrast(:,:,j);
                        else
                            s.minContrast(:,:,j) = [pad min(prev.minContrast(1:vStep:end-vStep,1+rel:end,j1),prev.minContrast(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                            s.maxContrast(:,:,j) = [pad max(prev.maxContrast(1:vStep:end-vStep,1+rel:end,j1),prev.maxContrast(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                        end
                        
                        s.minHalf(:,:,j) = [pad min(prev.diffs(1:vStep:end-vStep,1+rel:end,j1),prev.diffs(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                        s.maxHalf(:,:,j) = [pad max(prev.diffs(1:vStep:end-vStep,1+rel:end,j1),prev.diffs(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                        
                        if length > 1
                            s.minFourth(:,:,j) = [pad min(prev.minHalf(1:vStep:end-vStep,1+rel:end,j1),prev.minHalf(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                            s.maxFourth(:,:,j) = [pad max(prev.maxHalf(1:vStep:end-vStep,1+rel:end,j1),prev.maxHalf(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                        end
                        
                        if length > 2
                            s.minEigth(:,:,j) = [pad min(prev.minFourth(1:vStep:end-vStep,1+rel:end,j1),prev.minFourth(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                            s.maxEigth(:,:,j) = [pad max(prev.maxFourth(1:vStep:end-vStep,1+rel:end,j1),prev.maxFourth(1+vStep:vStep:end,1+rel-rel2:end-rel2,j2))];
                        end
                        
                        % increment angles
                        if j1 == j2 && i>1
                            j2 = j2+1;
                        else
                            j1 = j2;
                        end
                    % angle to the right, padding from the right boundary
                    else
                        s.sums(:,:,j) = [prev.sums(1:vStep:end-vStep,1:end-rel,j1)+prev.sums(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2) pad];
                        s.squares(:,:,j) = [prev.squares(1:vStep:end-vStep,1:end-rel,j1)+prev.squares(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2) pad];
                        s.diffs(:,:,j) = [prev.diffs(1:vStep:end-vStep,1:end-rel,j1)+prev.diffs(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2) pad];
                        s.diffSquares(:,:,j) = [prev.diffSquares(1:vStep:end-vStep,1:end-rel,j1)+prev.diffSquares(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2) pad];

                        if length <= prm.minContrast
                            s.minContrast(:,:,j) = [prev.diffs(1:vStep:end-vStep,1:end-rel,j1)+prev.diffs(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2) pad];
                            s.maxContrast(:,:,j) = s.minContrast(:,:,j);
                        else
                            s.minContrast(:,:,j) = [min(prev.minContrast(1:vStep:end-vStep,1:end-rel,j1),prev.minContrast(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                            s.maxContrast(:,:,j) = [max(prev.maxContrast(1:vStep:end-vStep,1:end-rel,j1),prev.maxContrast(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                        end
                        
                        s.minHalf(:,:,j) = [min(prev.diffs(1:vStep:end-vStep,1:end-rel,j1),prev.diffs(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                        s.maxHalf(:,:,j) = [max(prev.diffs(1:vStep:end-vStep,1:end-rel,j1),prev.diffs(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                        
                        if length > 1
                            s.minFourth(:,:,j) = [min(prev.minHalf(1:vStep:end-vStep,1:end-rel,j1),prev.minHalf(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                            s.maxFourth(:,:,j) = [max(prev.maxHalf(1:vStep:end-vStep,1:end-rel,j1),prev.maxHalf(1+vStep:vStep:end,1+rel2:end-rel+rel2,j2)) pad];
                        end
                        
                        % increment angles
                        if j1 == j2 && i>1
                            j1 = j1+1;
                        else
                            j2 = j1;
                        end
                    end
                end
                resp{i+1} = s;
            end
        end
        
        % Normalize the image values between 0 to 1.
        function I = normalize(I)
            I = I-min(I(:));
            I = I/max(I(:));
        end  
    end
end

