
% Compute the length of the major axis of a shape whose mask, centroid, and
% major axis orientation are given.  Assumes a concabe shape. Note that the
% MajorAxisLength property computed by regionprops is for the major axis of
% an approximating ellipse and not the actual axis of the shape.  
%
% mask - a mask of the nucleus containing all points in pointCA.
%
% centroid - a two element vetor whose first element is the
% x-coordinate of the nuceleus centroid and the second element is its
% y-coordinate.
%
% orientationDegrees - the angle formed by the major axis of the
% nucleus and the x-axis, measured from the x-axis to the major axis.


function [len x1 y1 x2 y2] = majorAxisLength(mask, centroid, orientationDegrees)
len = [];
x1 = [];
y1 = [];
x2 = [];
y2= [];

radiansPerDegree = pi / 180;

x = centroid(1);
y = centroid(2);

% Identify nucleus border pixels
erodedMask = imerode(mask, true(3));
border = mask & ~erodedMask;
[borderY, borderX] = find(border);

% Compute slope of line from centroid to each border pixel
borderSlope = (borderY - y) ./ (borderX - x);

axisSlope = tan(-orientationDegrees * radiansPerDegree);
diff = abs(borderSlope - axisSlope);

% diff contains nan values when inf is subtracted from inf.  Replace nan
%values with 0.
diff(isnan(diff)) = 0;
[~, order] = sort(diff);

% First end point
x1 = borderX(order(1));
y1 = borderY(order(1));

% Second end point
for i = 2:numel(order)
    x2 = borderX(order(i));
    y2 = borderY(order(i));
    % Make sure that second end point is at other end of the shape instead
    % of being next to first end point.  If the two points are on opposite
    % sides of the mask, then the distance between them is approximately
    % the sum of their distances to the centroid.  If the two points are on
    % the same side of the mask, then the distance from the centroid to the
    % further end point is approximately the sum of the distance from the
    % centroid to the closer point and the distance between the end points.
    ep1ep2Distance = distance(x1, y1, x2, y2);
    ep1centDistance = distance(x1, y1, x, y);
    ep2centDistance = distance(x2, y2, x, y);
    errOppSide = abs(ep1ep2Distance - (ep1centDistance + ep2centDistance));
    errSameSide = abs(max(ep1centDistance, ep2centDistance) - (min(ep1centDistance, ep2centDistance) + ep1ep2Distance));
%     fprintf('x1=%d y1==%d slope1=%f x2=%d y2=%d slope2=%f errOppSide=%f   errSameSide=%f\n', x1, y1, borderSlope(order(1)), x2, y2, borderSlope(order(i)), errOppSide, errSameSide);
    if errOppSide < errSameSide
        len = ep1ep2Distance;
        break;
    end
end
end


function d = distance(x1, y1, x2, y2)
d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

