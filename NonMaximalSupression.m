function [edgesFinal edgesImage] = NonMaximalSupression(edges,I,corrT)
    [m,n] = size(edges);
    [M,N] = size(I);
    
    edgesImage = zeros(M,N);
    binImage = zeros(M,N);
    edgesFinal = [];

    for i=1:m
        [curImage curThinImage] = getLineImage(edges(i,1),edges(i,2),edges(i,3),edges(i,4),I);
        L = edges(i,6);
        
        inter = curImage & binImage;
        corr = sum(inter(:))/(L);

        if corr < corrT && sum(curImage(:)) > 0
            edgesFinal = [edgesFinal ; edges(i,:) ];
            edgesImage = max(edgesImage,curThinImage*edges(i,5));
            binImage = binImage + curImage;
        end
    end
end

function [im imThin] = getLineImage(x0,y0,x1,y1,I)
    im = zeros(size(I));
    imThin = zeros(size(I));

    [M,N] = size(I);

    if(x0<=1 || y0<=1 || x1<=1 || y1<=1 || x0>=M || x1>=M || y0>=N || y1>=N)
        return;
    end
    
    dx = abs(x1-x0);
    dy = abs(y1-y0);
    
    shiftX = 1;
    shiftY = 1;
    
    if dx>dy
        shiftX = 0;
    else
        shiftY = 0;
    end

    if x0 < x1 
        sx = 1;
    else
        sx = -1;
    end

    if y0 < y1
        sy = 1;
    else
        sy = -1;
    end

    err = dx-dy;

    while 1>0
        imThin(x0,y0) = 1;
        im(x0,y0) = 1;
        im(x0+shiftX,y0+shiftY) = 1;
        im(x0-shiftX,y0-shiftY) = 1;
        %im(x0+2*shiftX,y0+2*shiftY) = 1;
        %im(x0+2*shiftX,y0+2*shiftY) = 1;


        if x0 == x1 && y0 == y1
            break;
        end

        e2 = 2*err;

        if e2 > -dy 
            err = err - dy;
            x0 = x0 + sx;
        end
        if e2 < dx
            err = err + dx;
            y0 = y0 + sy;
        end
    end
end
