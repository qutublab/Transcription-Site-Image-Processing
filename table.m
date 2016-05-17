function cluster3(fileName)
if nargin == 0
    fileName = 'out.csv';
end
[data varnames casenames] = tblread(fileName, ',');
varnamesCA = cell(size(varnames, 1), 1);
for i = 1:numel(varnamesCA)
   varnamesCA{i} = strtrim(varnames(i, :));
end


objectCountColumnName = 'Total Marker Object Count';
objectCountColumn = find(strcmp(varnamesCA, objectCountColumnName));
assert (numel(objectCountColumn) == 1, '%d occurrences of ''%s'' found', numel(objectCountColumn), objectCountColumnName);

radiusColumnName = 'Normalized Radius';
radiusColumn = find(strcmp(varnamesCA, radiusColumnName));
assert (numel(radiusColumn) == 1, '%d occurrences of ''%s'' found', numel(radiusColumn), radiusColumnName);

angleColumnName = 'Angle from Nucleus Major Axis';
angleColumn = find(strcmp(varnamesCA, angleColumnName));
assert (numel(angleColumn) == 1, '%d occurrences of ''%s'' found', numel(angleColumn), angleColumnName);

nucleusIndexColumnName = 'Nucleus Index';
nucleusIndexColumn = find(strcmp(varnamesCA, nucleusIndexColumnName));
assert (numel(nucleusIndexColumn) == 1, '%d occurrences of ''%s'' found', numel(nucleusIndexColumn), nucleusIndexColumnName);


radius = data(:, radiusColumn);
angle = data(:, angleColumn);
objectCount = data(:, objectCountColumn);
nucleusIndex = data(:, nucleusIndexColumn);



% For each point add another point 180 degrees away
radius = [radius; radius];
angle2 = rem(angle + 180, 360);
angle = [angle; angle2];
objectCount = [objectCount; objectCount];
nucleusIndex = [nucleusIndex; nucleusIndex];
markerFileName = [casenames; casenames];
% Convert angles to radians
angle = angle * (pi / 180);

% Use only the data from nuclei with 4 locations
idx4 = find(objectCount == 4);
radius = radius(idx4);
angle = angle(idx4);
nucleusIndex = nucleusIndex(idx4);
markerFileName = markerFileName(idx4, :);

% Convert to Cartesian coordinates
X = radius .* cos(angle);
Y = radius .* sin(angle);
points = [X, Y];

% Assume that 4 locations in same nucleus are consecutively listed in file.
% Break points into groups of 4.
numNuclei = numel(idx4) / 4;
markerLocations = cell(numNuclei, 1);
p = 1;
for i = 1:numNuclei
   markerLocations{i} = points(p:(p+3), :);
   p = p + 4;
%    figure, plot(X(p:(p+3)), Y(p:(p+3)), 's');
%    xlim([-1, 1])
%    ylim([-1, 1])
%    axis equal
end

distanceTable = zeros(numNuclei);
for i = 1:numNuclei
    loc1 = markerLocations{i};
    for j = 1:numNuclei
        loc2 = markerLocations{j};
        [~, D12] = knnsearch(loc1, loc2);
        [~, D21] = knnsearch(loc2, loc1);
        distanceTable(i, j) = sum(D12) + sum(D21);
    end
end

samePoints = distanceTable == 0;
distanceTable(samePoints) = Inf;
minDist = min(distanceTable(:));
closestPoints = distanceTable == minDist;
[R C] = find(closestPoints);
fprintf('Most similar nucleus pair(s) (distance=%f)\n', minDist);
for i = 1:numel(R)
    if R(i) < C(i)
        fprintf('nucleus %d in file %s and nucleus %d in file %s\n', ...
            nucleusIndex(R(i)), strtrim(markerFileName(R(i), :)), ...
            nucleusIndex(C(i)), strtrim(markerFileName(C(i), :)));
    end
end

outFileName = 'out-knn.csv';
fid = fopen(outFileName, 'w');
fprintf(fid, 'Nucleus');
for i = 1:numNuclei
    fprintf(fid, ',%d', i);
end
fprintf(fid, '\n');
for i = 1:numNuclei
    fprintf(fid, '%d', i);
    for j = 1:numNuclei
        fprintf(fid, ',%f', distanceTable(i, j));
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Wrote file %s\n', outFileName);
end

