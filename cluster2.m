function cluster2(fileName)

% radiusNormalization - 1: radius divided by distance of centroid to 
% boundary; 2: radius divided by nucleus major axis length
radiusNormalization = 1;

% Replicates parameter for kmeans function
kmeansReplicates = 1000;

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

% radiusColumnName = 'Normalized Radius';
% radiusColumn = find(strcmp(varnamesCA, radiusColumnName));
% assert (numel(radiusColumn) == 1, '%d occurrences of ''%s'' found', numel(radiusColumn), radiusColumnName);

angleColumnName = 'Marker Angle from Nucleus Major Axis';
angleColumn = find(strcmp(varnamesCA, angleColumnName));
assert (numel(angleColumn) == 1, '%d occurrences of ''%s'' found', numel(angleColumn), angleColumnName);

timeColumnName = 'Time in Minutes';
timeColumn = find(strcmp(varnamesCA, timeColumnName));
assert (numel(timeColumn) == 1, '%d occurrences of ''%s'' found', numel(timeColumn), timeColumnName);


% radius = data(:, radiusColumn);
[radius normalizationStr] = getMarkerDistance(varnamesCA, data, radiusNormalization);
angle = data(:, angleColumn);
objectCount = data(:, objectCountColumn);
timePoint = data(:, timeColumn);

for i = 1:max(objectCount)
    count = numel(find(objectCount == i));
    nucleiCount = count / i;
    fprintf('%d nuclei with %d locations; %d total locations\n', nucleiCount, i, count);
end
fprintf('\n');

numOriginalPoints = numel(radius);





% For each point add another point 180 degrees away
radius = [radius; radius];
angle2 = rem(angle + 180, 360);
angle = [angle; angle2];
% Convert angles to radians
angle = angle * (pi / 180);
objectCount = [objectCount; objectCount];
timePoint = [timePoint; timePoint];

% Since two points, one at 0 degrees one at 359 degrees are near each other
% switch to cartesian coordinates to refelect this for the Matlab kmeans
% function
X = radius .* cos(angle);
Y = radius .* sin(angle);
points = [X, Y];

marker = {{'p', 3}, {'o', 3}, {'^', 3}, {'s', 3}};

uniqueTimePoint = sort(unique(timePoint));
numUnique = numel(uniqueTimePoint);
timePointCount = zeros(numUnique, 1);
for i = 1:numUnique
   timePointCount(i) = sum(double(timePoint == uniqueTimePoint(i)));
end
legendLabels = cell(numUnique, 1);
legendMarkers = cell(numUnique, 1);
for i = 1:numUnique
    minutes = uniqueTimePoint(i);
    if minutes == 0
        legendLabels{i} = sprintf('Time not specified  (%d points)', timePointCount(i));
    else
        if rem(minutes, 60) == 0
            hrs = minutes / 60;
            legendLabels{i} = sprintf('%d hr  (%d points)', hrs, timePointCount(i));
        else
            legendLabels{i} = sprintf('%d min (%d points)', minutes, timePointCount(i));
        end
    end
    legendMarkers{i} = 's';
end
cmap = jet(numUnique);
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;

for i = 1:4
    markerSpec = marker{i};
    for j = 1:numUnique
        indices = (objectCount == i) & (timePoint == uniqueTimePoint(j));
        plot(points(indices, 1), points(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
    end
end
title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([-1, 1]);
    ylim([-1, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([-extremeX, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([-extremeY, extremeY]);
end
axis equal





numPoints = numel(radius);
maxK = round(sqrt(numPoints));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints, maxK);

bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
bestDistanceTable = [];
for k = 2:maxK
    [idx C sumd D] = kmeans(points, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(points, idx));
    fprintf('k: %d   mean silhoutte: %f\n', k, s);
    if s > bestS
        bestS = s;
        bestK = k;
        bestIdx = idx;
        bestCentroids = C;
        bestDistanceTable = D;
    end
end

fprintf('Best k between %d and %d: %d (mean silhouette: %f)\n', 2, maxK, bestK, bestS);



clusterSize = zeros(bestK, 1);
for i = 1:bestK
    clusterSize(i) = sum(double(bestIdx == i));
end


% figure, gscatter(points(:, 1), points(:, 2), bestIdx);
% axis square

cmap = jet(bestK);

legendLabels = cell(bestK, 1);
legendMarkers = cell(bestK, 1);
for i = 1:bestK
    legendLabels{i} = sprintf('Cluster %d  (%d points)', i, clusterSize(i));
    legendMarkers{i} = 's';
end
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;

for i = 1:4
    markerSpec = marker{i};
    for j = 1:bestK
        indices = (objectCount == i) & (bestIdx == j);
        plot(points(indices, 1), points(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
    end
end

title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; 'Clustered According to Cartesian Coordinates'; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([-1, 1]);
    ylim([-1, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([-extremeX, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([-extremeY, extremeY]);
end
axis equal







% Cluster only those locations that occur in groups of 4 within a nucleus
group4Indices = objectCount == 4;
points4 = points(group4Indices, :);
timePoint4 = timePoint(group4Indices);

uniqueTimePoint4 = sort(unique(timePoint4));

numUnique4 = numel(uniqueTimePoint4);
groupSize = zeros(numUnique4, 1);
for i = 1:numUnique4
    groupSize(i) = sum(double(timePoint4 == uniqueTimePoint4(i)));
end
legendLabels = cell(numUnique4, 1);
legendMarkers = cell(numUnique4, 1);
for i = 1:numUnique4
    minutes = uniqueTimePoint4(i);
    if minutes == 0
        legendLabels{i} = sprintf('Time not specified  (%d points)', groupSize(i));
    else
        if rem(minutes, 60) == 0
            hrs = minutes / 60;
            legendLabels{i} = sprintf('%d hr  (%d points)', hrs, groupSize(i));
        else
            legendLabels{i} = sprintf('%d min  (%d points)', minutes, groupSize(i));
        end
    end
    legendMarkers{i} = 's';
end
cmap = jet(numUnique4);
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;
markerSpec = marker{4};
for i = 1:numUnique4
    idx = (timePoint4 == uniqueTimePoint4(i));
    plot(points4(idx, 1), points4(idx, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
end
title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; '4 Transcription Sites per Nucleus'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([-1, 1]);
    ylim([-1, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([-extremeX, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([-extremeY, extremeY]);
end
axis equal




numPoints4 = size(points4, 1);
maxK = round(sqrt(numPoints4));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints4, maxK);

bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
bestDistanceTable = [];
for k = 2:maxK
    [idx C sumd D] = kmeans(points4, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(points4, idx));
    fprintf('k: %d   mean silhoutte: %f\n', k, s);
    if s > bestS
        bestS = s;
        bestK = k;
        bestIdx = idx;
        bestCentroids = C;
        bestDistanceTable = D;
    end
end
fprintf('Best k between %d and %d: %d (mean silhouette: %f)\n', 2, maxK, bestK, bestS);


clusterSize = zeros(bestK, 1);
for i = 1:bestK
   clusterSize(i) = sum(double(bestIdx == i));
end


cmap = jet(bestK);

markerSpec = marker{4};

legendLabels = cell(bestK, 1);
legendMarkers = cell(bestK, 1);
for i = 1:bestK
    legendLabels{i} = sprintf('Cluster %d  (%d points)', i, clusterSize(i));
    legendMarkers{i} = markerSpec{1};
end
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;

% newFig = true;
% for j = 1:bestK
%     indices = (bestIdx == j);
%     if newFig
%         figure, plot(points4(indices, 1), points4(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
%         hold on
%         newFig = false;
%     else
%         plot(points4(indices, 1), points4(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
%     end
% end

for j = 1:bestK
    indices = (bestIdx == j);
    plot(points4(indices, 1), points4(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
end


title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; 'Clustered According to Cartesian Coodinates'; '4 Transcription Sites per Nucleus'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([-1, 1]);
    ylim([-1, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([-extremeX, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([-extremeY, extremeY]);
end
axis equal
end
