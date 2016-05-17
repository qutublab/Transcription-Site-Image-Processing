
function nucleusMask = segmentDapi(fileName)
I = imread(fileName);
I = mat2gray(I);

% Segment by detecting edges
Gmag = imgradient(I);
M = im2bw(Gmag, graythresh(Gmag));        

% The next step will perform an imclose on individual objects.  Removing small
% objects especially if there a many of them keeps the processing time
% reasonable.
M = bwareaopen(M, 200);

% Close gaps individually so that two nuclei are not joined
[L numLabels] = bwlabel(M);
Lnz = L ~= 0;
acc = false(size(M));
se = strel('disk', 3, 0);
for k = 1:numLabels
    P = L == k;
%     P = imclose(P, se);
    % Remove pixels that overlap or abut another nucleus
    overlap = P & imdilate(acc | (L ~= k & Lnz), true(3));
    P = P & ~overlap;
    acc = acc | P;
end
        

M2 = imfill(acc, 'holes');     
M2 = imclearborder(M2);
M2 = bwareaopen(M2, 2000);
        
% Keep only non-clusters
[L numLabels] = bwlabel(M2);
props = regionprops(L, 'Solidity');
nucleusMask = false(size(L));
for i = 1:numLabels
    if props(i).Solidity >= 0.95
        nucleusMask(L == i) = true;
    end
end

% r = I;
% g = I;
% b = I;
% brdr = nucleusMask & ~imerode(nucleusMask, true(3));
% r(brdr) = 0; g(brdr) = 0; b(brdr) = 1;
% figure, imshow(cat(3, r, g, b));
% title(fileName);

end
