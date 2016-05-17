function BW = drawline2(x, y, x2, y2, BW)
[numRows, numCols] = size(BW);

x = round(x);
y = round(y);
x2 = round(x2);
y2 = round(y2);
slope = (y2 - y) / (x2 - x);
if abs(slope) <= 1
    % Use y = mx + b
    yIntercept = y - (x * slope);
    xRange = x:sign(x2-x):x2;
    yRange = round((xRange * slope) + yIntercept);
else
    % Use x = my + b
    slope = 1 / slope;
    xIntercept = x - (y * slope);
    yRange = y:sign(y2-y):y2;
    xRange = round((yRange * slope) + xIntercept);
end

linearIndices = sub2ind([numRows numCols], yRange, xRange);
BW(linearIndices) = true;



end