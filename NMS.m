function EdgeListFinal = NMS(EdgeList,LeechMask,sigma,Frac,image)
    % Edge list table
    % 1. level number
    % 2. "1" for a vertical edge "3" for horizontal edge
    % 3. x1 first x coordinate 
    % 4. y1 first y coordinate
    % 5. x2 second x coordinate
    % 6. y2 second y coordinate
    % 7. orientation number

    EndLevel = max(EdgeList(:,1));
    StartLevel = min(EdgeList(:,1));

    EdgeListFinal = [EdgeList zeros(size(EdgeList,1),5)];

    %ResponsePixelMap = zeros(size(image)); 
    Raster = zeros(size(image));
    SoftResponsePixelMap = zeros(size(image));
    %LeechMap = zeros(size(image));
    %ResponsePixelMap_NoLeech = zeros(size(image));
    %LevelMap = zeros(size(image));

    %
    [len1 len2] = size(image);
    logN = log(len1 * len2);
    sumabsLeechMask = sum(abs(LeechMask));
    % Prepare Responses obtained by LeechMask in advance
    ResponseVertical_LeechMask = conv2(double(image),LeechMask,'same'); 
    ResponseHorizontal_LeechMask = conv2(double(image),LeechMask','same');


    len = 2.^(0:1:EndLevel-1);
    ThreshValue = sigma*sqrt(2 * logN ./ (sumabsLeechMask .* len));

    for Level=EndLevel:-1:StartLevel
        % remark: the number of directions per length is 2^Level+1 (where Level = log_2(length)+1) 
        % Since we have two sets per a certain length, i.e. shift and
        % non-shift, and we would like to scan each of them seprately, we
        % enumarte them by ...
        num_of_directions = 2*(2^Level + 1); 
        % factor that moves L2 norm to Linf norm.
        % We don't have to move from L2 to Linf,
        % the calculation of the tail length is already Linf
        % directions1 = (1:1:num_of_directions/2);
        % directions = [directions1 directions1];
        % correct_tail_length = len(Level) ./ sqrt((directions - len(Level)-1).^2 + len(Level)^2);
        %%%%%%% VERTICAL EDGES %%%%%%%
        % scan the vertical edges at each level according to their direction 
        for t=1:num_of_directions
            IndexVertical = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 1) & (EdgeList(:,7) == t));
            num_of_pairs = size(IndexVertical,1);
            if (num_of_pairs == 0)
                continue;
            end

            x1 = EdgeList(IndexVertical,3);
            y1 = EdgeList(IndexVertical,4);
            x2 = EdgeList(IndexVertical,5);
            y2 = EdgeList(IndexVertical,6);
    %
            color = (1:1:num_of_pairs);
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','Custom');%setting border color to custom.
            set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
            pts=[x1, y1, x2, y2];
            PixelMap = zeros(size(image));
            PixelMap = step(h,PixelMap,pts);%performs the drawing.
            % Calculate edge length at once 
            count_edge_length = hist(PixelMap(:),(0:1:num_of_pairs));
            % Calculate intersections with previous Raster at once
            joint_Raster = Raster .* PixelMap;
            count_intersect = hist(joint_Raster(:),(0:1:num_of_pairs));
            % Calculate non-intersecting parts (tails) 
            non_joint_Raster = ~Raster .* PixelMap;
            tail_length = hist(non_joint_Raster(:),(0:1:num_of_pairs));
            tail_thresh =  sigma*sqrt(2 * logN ./ (sumabsLeechMask * tail_length));

    %
            for I=1:num_of_pairs
    % 
               if (count_intersect(I+1) / count_edge_length(I+1) <= Frac)  
                    % Leech test
                    non_joint_Raster_TMP = (non_joint_Raster == I);
    %                tail_mask = conv2(double(non_joint_Raster_TMP),LeechMask,'same');
    %                tail_response = abs(sum(sum(tail_mask.*image)) / (sumabsLeechMask * tail_length(I+1)));
                    tail_response = abs(sum(sum(non_joint_Raster_TMP.*ResponseVertical_LeechMask))) / ... 
                        (sumabsLeechMask * tail_length(I+1));
                    if (tail_response < tail_thresh(I+1))
                         EdgeListFinal(IndexVertical(I),end) = -1; %Leech
                    else
                        EdgeListFinal(IndexVertical(I),end-4) = len(Level);
                        EdgeListFinal(IndexVertical(I),end-3) = tail_length(I+1);
                        EdgeListFinal(IndexVertical(I),end-2) = tail_response;
                        EdgeListFinal(IndexVertical(I),end-1) = tail_thresh(I+1);
                        EdgeListFinal(IndexVertical(I),end) = 1; %Edge
                    end
                end
            end
            % Accumulate Soft Map of Responses
            IndexVerticalEdges = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 1) & (EdgeList(:,7) == t) & (EdgeListFinal(:,end) == 1));
            if ~isempty(IndexVerticalEdges)
                x1 = EdgeList(IndexVerticalEdges,3);
                y1 = EdgeList(IndexVerticalEdges,4);
                x2 = EdgeList(IndexVerticalEdges,5);
                y2 = EdgeList(IndexVerticalEdges,6);    
                % The response is the color  
                color = EdgeList(IndexVerticalEdges,8) / ThreshValue(Level);
                h = vision.ShapeInserter;
                release(h);%allows changing ShapeInserter object properties.
                set(h,'Shape','Lines');
                set(h,'BorderColor','Custom');%setting border color to custom.
                set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
                pts=[x1, y1, x2, y2];
                SoftResponsePixelMap_TMP = zeros(size(image));
                SoftResponsePixelMap_TMP = step(h,SoftResponsePixelMap_TMP,pts);%performs the drawing.
                SoftResponsePixelMap = max(double(SoftResponsePixelMap),double(SoftResponsePixelMap_TMP));
            end
        end
        % Update the Raster only after scanning the whole set of responses of a
        % certain level, by convolution
        % Use Vision Toolbox to prepare RasterTMP, Cheaper !!!    
        IndexVertical = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 1));
        num_of_pairs = size(IndexVertical,1);
        if (num_of_pairs > 0)
            x1 = EdgeList(IndexVertical,3);
            y1 = EdgeList(IndexVertical,4);
            x2 = EdgeList(IndexVertical,5);
            y2 = EdgeList(IndexVertical,6);    
            %
            RasterTMP = zeros(size(image));
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','White');%draws white color on black image.
            pts=[x1, y1, x2, y2];
            RasterTMP = step(h,RasterTMP,pts);%performs the drawing.
            Raster = Raster | conv2(double(RasterTMP),ones(1,3),'same');   
        end


    % scan the horizontal edges at each level according to their direction     
         for t=1:num_of_directions
             %%%%%%% HORIZONTAL EDGES %%%%%%%
             IndexHorizontal = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 3) & (EdgeList(:,7) == t));
             num_of_pairs = size(IndexHorizontal,1);
             if (num_of_pairs == 0)
                continue;
             end
    %
             x1 = EdgeList(IndexHorizontal,3);
             y1 = EdgeList(IndexHorizontal,4);
             x2 = EdgeList(IndexHorizontal,5);
             y2 = EdgeList(IndexHorizontal,6);
    %
            color = (1:1:num_of_pairs);
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','Custom');%setting border color to custom.
            set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
            pts=[x1, y1, x2, y2];
            PixelMap = zeros(size(image));
            PixelMap = step(h,PixelMap,pts);%performs the drawing.
            % Calculate edge length at once 
            count_edge_length = hist(PixelMap(:),(0:1:num_of_pairs));
            % Calculate intersections with previous Raster at once
            joint_Raster = Raster .* PixelMap;
            count_intersect = hist(joint_Raster(:),(0:1:num_of_pairs));
            % Calculate non-intersecting parts (tails) 
            non_joint_Raster = ~Raster .* PixelMap;
            tail_length = hist(non_joint_Raster(:),(0:1:num_of_pairs));
            tail_thresh =  sigma*sqrt(2 * logN ./ (sumabsLeechMask * tail_length));
    %
            for I=1:num_of_pairs
                  if (count_intersect(I+1) / count_edge_length(I+1) <= Frac)           
                    % Leech test
                    non_joint_Raster_TMP = (non_joint_Raster == I);
    %                tail_mask = conv2(double(non_joint_Raster_TMP),LeechMask','same');
    %                tail_response = abs(sum(sum(tail_mask.*image)) / (sumabsLeechMask * tail_length(I+1)));
                    tail_response = abs(sum(sum(non_joint_Raster_TMP.*ResponseHorizontal_LeechMask))) / ... 
                        (sumabsLeechMask * tail_length(I+1));
                    if (tail_response < tail_thresh(I+1))
                        EdgeListFinal(IndexHorizontal(I),end) = -1; % Leech
                    else
                        EdgeListFinal(IndexHorizontal(I),end-4) = len(Level);
                        EdgeListFinal(IndexHorizontal(I),end-3) = tail_length(I+1);
                        EdgeListFinal(IndexHorizontal(I),end-2) = tail_response;
                        EdgeListFinal(IndexHorizontal(I),end-1) = tail_thresh(I+1);
                        EdgeListFinal(IndexHorizontal(I),end) = 1; %Edge
                    end
                 end
            end
            % Accumulate Soft Map of Responses
            IndexHorizontalEdges = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 3) & (EdgeList(:,7) == t) & (EdgeListFinal(:,end) == 1));
            if ~isempty(IndexHorizontalEdges)
                x1 = EdgeList(IndexHorizontalEdges,3);
                y1 = EdgeList(IndexHorizontalEdges,4);
                x2 = EdgeList(IndexHorizontalEdges,5);
                y2 = EdgeList(IndexHorizontalEdges,6);    
                % The response is the color  
                color = EdgeList(IndexHorizontalEdges,8) / ThreshValue(Level);
                h = vision.ShapeInserter;
                release(h);%allows changing ShapeInserter object properties.
                set(h,'Shape','Lines');
                set(h,'BorderColor','Custom');%setting border color to custom.
                set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
                pts=[x1, y1, x2, y2];
                SoftResponsePixelMap_TMP = zeros(size(image));
                SoftResponsePixelMap_TMP = step(h,SoftResponsePixelMap_TMP,pts);%performs the drawing.
                SoftResponsePixelMap = max(double(SoftResponsePixelMap),double(SoftResponsePixelMap_TMP));
            end
         end

    % Update the Raster only after scanning the whole set of responses of a
    % certain level, by convolution
    % Use Vision Toolbox to prepare RasterTMP, Cheaper !!!
        IndexHorizontal = find((EdgeList(:,1) == Level) & (EdgeList(:,2) == 3));
        num_of_pairs = size(IndexHorizontal,1);
        if (num_of_pairs > 0)
            x1 = EdgeList(IndexHorizontal,3);
            y1 = EdgeList(IndexHorizontal,4);
            x2 = EdgeList(IndexHorizontal,5);
            y2 = EdgeList(IndexHorizontal,6);    
            RasterTMP = zeros(size(image));
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','White');%draws white color on black image.
            pts=[x1, y1, x2, y2];
            RasterTMP = step(h,RasterTMP,pts);%performs the drawing.
            Raster = Raster | conv2(double(RasterTMP),ones(3,1),'same'); 
        end 

    end
    % Generate Binary Map of Edges Before Leech Removal
    ResponsePixelMap = zeros(size(image));
    IndexPotentialEdges = find(EdgeListFinal(:,end) ~= 0); 
    if ~isempty(IndexPotentialEdges)
            x1 = EdgeList(IndexPotentialEdges,3);
            y1 = EdgeList(IndexPotentialEdges,4);
            x2 = EdgeList(IndexPotentialEdges,5);
            y2 = EdgeList(IndexPotentialEdges,6);
    %
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','White');%draws white color on black image.
            pts=[x1, y1, x2, y2];
            ResponsePixelMap = zeros(size(image));
            ResponsePixelMap = step(h,ResponsePixelMap,pts);%performs the drawing.
    end
    % Generate Binary Map of Leeches
    LeechMap = zeros(size(image));
    IndexLeech = find(EdgeListFinal(:,end) == -1); 
    if ~isempty(IndexLeech)
            x1 = EdgeList(IndexLeech,3);
            y1 = EdgeList(IndexLeech,4);
            x2 = EdgeList(IndexLeech,5);
            y2 = EdgeList(IndexLeech,6);
    %
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','White');%draws white color on black image.
            pts=[x1, y1, x2, y2];
            LeechMap = zeros(size(image));
            LeechMap = step(h,LeechMap,pts);%performs the drawing.
    end
    % Generate Binary Map of Edges, After Leech Removal
    ResponsePixelMap_NoLeech = zeros(size(image));
    IndexEdges = find(EdgeListFinal(:,end) == 1); 
    if ~isempty(IndexEdges)
            x1 = EdgeList(IndexEdges,3);
            y1 = EdgeList(IndexEdges,4);
            x2 = EdgeList(IndexEdges,5);
            y2 = EdgeList(IndexEdges,6);
    %
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','White');%draws white color on black image.
            pts=[x1, y1, x2, y2];
            ResponsePixelMap_NoLeech = zeros(size(image));
            ResponsePixelMap_NoLeech = step(h,ResponsePixelMap_NoLeech,pts);%performs the drawing.
    end
    % Generate Level Map, colors according to the level
    LevelMap = zeros(size(image));
    if ~isempty(IndexEdges)
            color = (EdgeList(IndexEdges,1));
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','Custom');%setting border color to custom.
            set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
            pts=[x1, y1, x2, y2];
            LevelMap = zeros(size(image));
            LevelMap = step(h,LevelMap,pts);%performs the drawing.
    end

    if 1
        if 1
    %       figure; imagesc(image); axis('image'); colormap(gray); title('original image');
            ResponsePixelMap_overlay=MarkBib2im(image,ResponsePixelMap);
            figure; imagesc(ResponsePixelMap_overlay); axis('image'); title('Response Map BEFORE Leech Removal AFTER inter-level suppression');
            LeechMap_overlay = MarkBib2im(image,LeechMap);
            figure; imagesc(LeechMap_overlay); axis('image'); title('Leech Map');
    %       ResponsePixelMap_overlay = MarkBib2im(image,~LeechMap .* ResponsePixelMap);
    %       figure; imagesc(ResponsePixelMap_overlay); axis('image'); title('Response Map AFTER Leech Removal');
        end
        ResponsePixelMap_overlay = MarkBib2im(image,ResponsePixelMap_NoLeech);
        if 1
    %        figure; imagesc(ResponsePixelMap_overlay); axis('image'); title('Response Map AFTER NMS and Leech Removal');
            figure; imagesc(ResponsePixelMap_overlay); axis('image'); title('Response Map','FontSize',20);
            figure; imagesc(SoftResponsePixelMap); axis('image'); colormap(gray); title('Soft Response','FontSize',20);
        end
        IndexedLevelMap = LevelMap;
        if 1
            figure; imagesc(LevelMap); axis('image'); title('Level Map','FontSize',20);
    %        storecolormap = [0 0.2 0.4 0.6 0.7 0.8 0.9 0.95 1; 0 0.2 0.3 0.4, 0.4, 0.3, 0.2, 0 0; 1 0.95 0.9 0.8 0.7 0.6 0.4 0.2 0]';
    %        colormap(storecolormap);
            storecolormap = colormap;
            IndexedLevelMap = ind2rgb(LevelMap*10,storecolormap); % general
    % for Protrack
    %        load EdgesColorMap; % use colormapeditor to edit the color map and save to EdgesColorMap
    %        colormap(EdgesColorMap);
        end
    end

    % profile viewer;
    % These are the outputs !!
    EdgeListFinal = EdgeListFinal(EdgeListFinal(:,end) > 0,:);
    EdgeListFinal = EdgeListFinal(:,1:end-1);
    ee = SoftResponsePixelMap; % Soft Map
    eeb = ResponsePixelMap_NoLeech; % Binary Map
    eeb_overlay = ResponsePixelMap_overlay; % Binary Map overlaid on original image
    %IndexedLevelMap=ind2rgb(LevelMap*10,EdgesColorMap); % Indexed LevelMap "HexagonResults1"
    %IndexedLevelMap=ind2rgb(LevelMap*9,EdgesColorMap);   % "HexagonResults2"
    %IndexedLevelMap = ind2rgb(LevelMap*10,storecolormap); % general
    EdgeTable = EdgeListFinal;
    %% Draw two images, one for the SIGNED horizontal responses and the other for the SIGNED vertical responses 
    % prepare two maps which will be used for fiber detection 
    VerticalResponse = zeros(size(image));
    IndexEdges = find(EdgeListFinal(:,2) == 1); 
    if ~isempty(IndexEdges)
    %
            x1 = EdgeListFinal(IndexEdges,3);
            y1 = EdgeListFinal(IndexEdges,4);
            x2 = EdgeListFinal(IndexEdges,5);
            y2 = EdgeListFinal(IndexEdges,6);
    %        
            color = (EdgeListFinal(IndexEdges,9));
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','Custom');%setting border color to custom.
            set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
            pts=[x1, y1, x2, y2];
            VerticalResponse = step(h,VerticalResponse,pts);%performs the drawing.
    end
    %

    HorizontalResponse = zeros(size(image));
    IndexEdges = find(EdgeListFinal(:,2) == 3); 
    if ~isempty(IndexEdges)
    %
            x1 = EdgeListFinal(IndexEdges,3);
            y1 = EdgeListFinal(IndexEdges,4);
            x2 = EdgeListFinal(IndexEdges,5);
            y2 = EdgeListFinal(IndexEdges,6);
    %        
            color = (EdgeListFinal(IndexEdges,9));
            h = vision.ShapeInserter;
            release(h);%allows changing ShapeInserter object properties.
            set(h,'Shape','Lines');
            set(h,'BorderColor','Custom');%setting border color to custom.
            set(h,'CustomBorderColor',color);%sets line color,color is an R element vector.
            pts=[x1, y1, x2, y2];
            HorizontalResponse = step(h,HorizontalResponse,pts);%performs the drawing.
    end
end