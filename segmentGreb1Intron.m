

function [intronMask preMask] = segmentGreb1Intron(fileName)

J = imread(fileName);
J = mat2gray(J);
preMask = im2bw(J, graythresh(J));
intronMask = imclose(preMask, strel('disk', 7, 0));
intronMask = bwareaopen(intronMask, 10);
 
end
