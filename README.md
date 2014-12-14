FaintStraightEdgeDetection
==========================

Detector of faint straight edges in noisy images. 
It runs fast, but do not detect every edge. 
If you want to invest more run time and to get better quality of edge detection I suggest you to try my triangles solutions:
https://github.com/NatiOfir/TrianglesEdgeDetection.git

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm
