FaintStraightEdgeDetection
==========================

Detector of faint straight edges in noisy images. This work is described in the TPAMI 2020 paper:
https://ieeexplore.ieee.org/document/8607091


It runs fast, but do not detect every edge.
If you want to invest more run time and to get better quality of edge detection I suggest you to try my triangle solutions:
https://github.com/NatiOfir/TrianglesEdgeDetection.git

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm

Good luck!
