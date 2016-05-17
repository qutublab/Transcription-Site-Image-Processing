function BW = drawline(x, y, slope, BW)
[numRows, numCols] = size(BW);

x = round(x);
y = round(y);
if abs(slope) <= 1
    % Use y = mx + b
    yIntercept = y - (x * slope);
    xRange = 1:numCols;
    yRange = round((xRange * slope) + yIntercept);
    inMask = (yRange >= 1) & (yRange <= numRows);
else
    % Use x = my + b
    slope = 1 / slope;
    xIntercept = x - (y * slope);
    yRange = 1:numRows;
    xRange = round((yRange * slope) + xIntercept);
    inMask = (xRange >= 1) & (xRange <= numCols);
end
xRange = xRange(inMask);
yRange = yRange(inMask);

linearIndices = sub2ind([numRows numCols], yRange, xRange);
BW(linearIndices) = true;



end