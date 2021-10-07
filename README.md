FaintStraightEdgeDetection
==========================

Detector of faint straight edges in noisy images. This work is described in the TPAMI 2020 paper "On Detection of Faint Edges in Noisy Images":
https://ieeexplore.ieee.org/document/8607091


It runs fast but detects only straight edges.
If you want to invest more run time and to get better quality of edge detection I suggest you to try the triangle solutions for example:
https://github.com/NatiOfir/TrianglesEdgeDetection.git

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm

Good luck!
