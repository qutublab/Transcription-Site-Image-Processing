

% Compute the polar coordinates of points within a nucleus from their
% cartesian coordinates in the cell array pointCA.  
%
% pointCA - a cell array of two-element vectors where the first element is
% the x-coordinate of a point and the second element is the y-ccordinate of
% the point.
%
% nucleusMask - a mask of the nucleus containing all points in pointCA.
%
% nucleusCentroid - a two element vetor whose first element is the
% x-coordinate of the nuceleus centroid and the second element is its
% y-coordinate.
%
% nucleusOrientationDegrees - the angle formed by the major axis of the
% nucleus and the x-axis, measured from the x-axis to the major axis.
%
% rThetaArr - an array with n rows and 3 columns where n is the number of
% points of pointCA.  Each row contains the distance from the nucleus
% centroid to the point, the angle (degrees) from the major axis of the
% nucleus to the point (positive angles indicate a counter-clockwise
% rotation), and the distance of a line segment from the nucleus centroid
% passing through the point and terminating at the nucleus border.


function rThetaArr = polar(pointCA, nucleusMask, nucleusCentroid, nucleusOrientationDegrees)
degreesPerRadian = 180 / pi;

% Display results, part 1 of 5
% nuc = false(size(nucleusMask));
% pt = nuc;
% brdr = nuc;
% lines = nuc;
% maj = nuc;
% ref = nuc;

nucleusX = nucleusCentroid(1);
nucleusY = nucleusCentroid(2);

% Display results, part 2 of 5
% nuc(round(nucleusY), round(nucleusX)) = true;
% maj = drawline(nucleusX, nucleusY, tan(-nucleusOrientationDegrees/degreesPerRadian), maj);
% refAngle = -nucleusOrientationDegrees + 90;
% ref = drawline(nucleusX, nucleusY, tan(refAngle/degreesPerRadian), ref);
% refAngle = -nucleusOrientationDegrees + 45;
% ref = drawline(nucleusX, nucleusY, tan(refAngle/degreesPerRadian), ref);
% refAngle = -nucleusOrientationDegrees -45;
% ref = drawline(nucleusX, nucleusY, tan(refAngle/degreesPerRadian), ref);




% Identify nucleus border pixels
erodedNucleusMask = imerode(nucleusMask, true(3));
nucleusBorder = nucleusMask & ~erodedNucleusMask;
[nucleusBorderY, nucleusBorderX] = find(nucleusBorder);

% Compute slope of line from nucleus centroid to each border pixel
borderSlope = (nucleusBorderY - nucleusY) ./ (nucleusBorderX - nucleusX);

numPoints = numel(pointCA);
rThetaArr = zeros(numPoints, 3);
for j = 1:numPoints
    % Find the nucleus border pixel that is on the ray that starts at the
    % nucleus centroid and passes through the given point.  This is
    % done by finding the nucleus border pixel with the closest slope
    % (relative to the nucleus centroid) to the slope of the given point.

    point = pointCA{j};
    x = point(1);
    y = point(2);

% Display results, part 3 of 5
% pt(round(y), round(x)) = true;

    deltaY = y - nucleusY;
    deltaX = x - nucleusX;
    slope = deltaY / deltaX;
    
    diff = abs(slope - borderSlope);
    % diff contains nan values when inf is subtracted from inf.  Replace
    % nan values with 0.
    diff(isnan(diff)) = 0;
    [~, order] = sort(diff);

    % Note that the border pixel with the best slope could be 180 degrees away
    % from the desired border pixel because of how the slope is calculated.
    % Thus we must examine pixel coordinates.  Compare three distances:
    % distA: the distance between the nucleus centroid and the given point;
    % distB: the distance between the given point and the
    % border pixel; and distC: the distance between the nucleus centroid
    % and the border pixel.  For the desired border pixel, the given point
    % is between it and the nucleus centroid, i.e.
    % distA + distB = distC. If the border pixel is 180 degrees from the
    % desired border pixel, then the nucleus centroid is between it and the
    % given point, i.e. distB = distA + distC.
    distA = distance(nucleusX, nucleusY, x, y);
    borderPixelFound = false;
    for k = 1:numel(order)
        borderX = nucleusBorderX(order(k));
        borderY = nucleusBorderY(order(k));
        distB = distance(x, y, borderX, borderY);
        distC = distance(nucleusX,  nucleusY, borderX,  borderY);
        % Determine if the pixel is in the desired location or is 180
        % degrees away from the desired location.
        err = abs(distC - (distA + distB));
        err180 = abs(distB - (distA + distC));
        if err < err180
            borderPixelFound = true;
            break;
        end
    end
    if ~borderPixelFound
        fprintf('Nucleus centroid: (%f, %f) Intron point: (%f, %f)\n', nucleusX, nucleusY, x, y);
        se = strel('disk', 5, 0);
        nucleusCenter = false(size(nucleusMask));
        nucleusCenter(round(nucleusY), round(nucleusX)) = true;
        nucleusCenter = imdilate(nucleusCenter, se);
        intronPoint = false(size(nucleusMask));
        intronPoint(round(y), round(x)) = true;
        intronPoint = imdilate(intronPoint, se);
        r = nucleusMask & nucleusCenter & ~intronPoint;
        g = nucleusMask & ~nucleusCenter & intronPoint;
        b = nucleusMask & ~nucleusCenter & ~intronPoint;
        figure, imshow(double(cat(3, r, g, b)));
        figure, imshow(nucleusMask);
        figure, imshow(nucleusCenter);
        figure, imshow(intronPoint);
    end
    assert(borderPixelFound, '[polar.m ] Border pixel not found for point: %d', j);

% Display results, part 4 of 5
% lines = drawline(nucleusX, nucleusY, slope, lines);
% brdr(borderY, borderX) = true;

    % Normalize polar radius
    %r = distA / distC;
    
    % Theta, the angle component of the polar coordinates of the given
    % point, is to be measured from the major axis of the nucleus.  Let
    % alpha be the angle from the x-axis to the major axis of the nucleus.
    % Let beta be the angle from the x-axis to the line formed by the
    % nucleus centroid to the given point.  Hence
    % beta = theta + alpha, and theta = beta - alpha
    alphaDegrees = -nucleusOrientationDegrees;
    betaTangent = slope;
    betaDegrees = atan(betaTangent) * degreesPerRadian;
    % Note that the range of atan is -pi/2 to pi/2.  So, if the point is in
    % quadrants II or III, betaDegrees must be adjusted
    if deltaX < 0
        if deltaY > 0      % Quadrant II
            betaDegrees = betaDegrees + 180;
        elseif deltaY < 0  % Quadrant III
            betaDegrees = betaDegrees - 180;
        end
    end
        
    thetaDegrees = betaDegrees - alphaDegrees;
    % Force theta to be between 0 and 360 degrees
    thetaDegrees = rem(thetaDegrees, 360);
    if thetaDegrees < 0;
        thetaDegrees = thetaDegrees + 360;
    end
%     fprintf('alpha=%f  beta=%f  theta=%f   r=%f\n', alphaDegrees, betaDegrees, thetaDegrees, r);
    
%     rThetaArr(j,:) = [r, thetaDegrees];
    rThetaArr(j,:) = [distA, distC, thetaDegrees];
end

% Display results, part 5 of 5
% nuc = imdilate(nuc, true(5));
% pt = imdilate(pt, true(5));
% brdr = imdilate(brdr, true(5));
% r = nucleusMask;
% g = nucleusMask;
% b = nucleusMask;
% r(ref) = true;
% g(ref) = true;
% b(ref) = false;
% r(lines) = false;
% g(lines) = false;
% b(lines) = false;
% r(maj) = true;
% g(maj) = false;
% b(maj) = false;
% 
% r(nuc) = true;
% g(nuc) = false;
% b(nuc) = false;
% r(pt) = false;
% g(pt) = true;
% b(pt) = false;
% r(brdr) = false;
% g(brdr) = false;
% b(brdr) = true;
% 
% rgb = double(cat(3, r(end:-1:1,:), g(end:-1:1,:), b(end:-1:1,:)));
% figure, imshow(rgb);


end


function d = distance(x1, y1, x2, y2)
d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

