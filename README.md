FaintStraightEdgeDetection
==========================

Detector of faint straight edges in noisy images.

This code implemants the algorithm that is described in the following paper:
M. Galun, R. Basri and A. Brandt
Multiscale edge detection and fiber enhancement using differences of oriented means 
IEEE International Conference on Computer Vision, Rio De Janeiro, ICCV-07 (2007)
http://www.wisdom.weizmann.ac.il/~/meirav/EdgesGalunBasriBrandt.pdf

It runs fast, but do not detect every edge.
If you want to invest more run time and to get better quality of edge detection I suggest you to try my triangle solutions:
https://github.com/NatiOfir/TrianglesEdgeDetection.git

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm
