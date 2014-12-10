function R = runIm( I,sigma )
        im = Image(I);        
        im = im.orientedMeans;
        im = im.detectEdges(0,sigma,sigma);
        im = im.nonMaximalSupression(true);
        R = im.edgesImage;
end

