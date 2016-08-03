function cluster4(fileName)

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


% Cluster only those locations that occur in groups of 3 within a nucleus
group3Indices = objectCount == 3;
points3 = points(group3Indices, :);
timePoint3 = timePoint(group3Indices);

uniqueTimePoint3 = sort(unique(timePoint3));
numUnique3 = numel(uniqueTimePoint3);
groupCount = zeros(numUnique3, 1);
for i = 1:numUnique3
   groupCount(i) = sum(double(timePoint3 == uniqueTimePoint3(i))); 
end

legendLabels = cell(numUnique3, 1);
legendMarkers = cell(numUnique3, 1);
for i = 1:numUnique3
    minutes = uniqueTimePoint3(i);
    if minutes == 0
        legendLabels{i} = sprintf('Time not specified  (%d points)', groupCount(i));
    else
        if rem(minutes, 60) == 0
            hrs = minutes / 60;
            legendLabels{i} = sprintf('%d hr  (%d points)', hrs, groupCount(i));
        else
            legendLabels{i} = sprintf('%d min (%d points)', minutes, groupCount(i));
        end
    end
    legendMarkers{i} = '^';
end
cmap = jet(numUnique3);
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;
markerSpec = marker{3};
for i = 1:numUnique3
    idx = (timePoint3 == uniqueTimePoint3(i));
    plot(points3(idx, 1), points3(idx, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
end
title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; '3 Transcription Sites per Nucleus'});
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




numPoints3 = size(points3, 1);
maxK = round(sqrt(numPoints3));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints3, maxK);


bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
bestDistanceTable = [];
for k = 2:maxK
    [idx C sumd D] = kmeans(points3, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(points3, idx));
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
   clusterSize(i) =  sum(double(bestIdx == i));
end

cmap = jet(bestK);

markerSpec = marker{3};

legendLabels = cell(bestK, 1);
legendMarkers = cell(bestK, 1);
for i = 1:bestK
    legendLabels{i} = sprintf('Cluster %d  (%d points)', i, clusterSize(i));
    legendMarkers{i} = markerSpec{1};
end
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;


for j = 1:bestK
    indices = (bestIdx == j);
    plot(points3(indices, 1), points3(indices, 2), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
end


title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; 'Clustered According to Cartesian Coordinates'; '3 Transcription Sites per Nucleus'});
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



% Plot cluster data as 4 concentric rings where borders represent the distance
% of the 4 farthest points in each quartile

numRings = 4;
cmap = jet(numRings);
h = figure;
for i = 1:bestK
    inCluster = bestIdx == i;
    clusterXY = points3(inCluster, :);
    clusterSize = size(clusterXY, 1);
    centroidX = bestCentroids(i, 1);
    centroidY = bestCentroids(i, 2);
    % Calculate distance of each point in cluster from cluster center
    distanceFromClusterCenter = sqrt((clusterXY(:, 1) - centroidX).^2 + (clusterXY(:, 2) - centroidY).^2);
    distanceFromClusterCenter = sort(distanceFromClusterCenter);
    for r = numRings:-1:1
        d = distanceFromClusterCenter(round(clusterSize * (r / numRings)));
        drawDisk(centroidX, centroidY, d, cmap(r, :));
    end
end
title({'Cartesian Coordinates of Intron Locations with 180 Degree Replicates'; normalizationStr; 'Clustered According to Cartesian Coordinates'; '3 Transcription Sites per Nucleus'; 'Colored regions within each cluster contain equal numbers of points'});
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
