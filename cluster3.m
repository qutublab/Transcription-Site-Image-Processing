
function cluster3(fileName)

% radiusNormalization - 1: radius divided by distance of centroid to 
% boundary; 2: radius divided by nucleus major axis length
radiusNormalization = 1;

% Replicates parameter for kmeans function
kmeansReplicates = 1; %1000;


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

% Rotate points by 90, 180, or 270 degrees so that the are in the first
% quadrant
angle = rem(angle, 90);

% Shapes and sizes of plot markers
marker = {{'p', 3}, {'o', 3}, {'^', 3}, {'s', 3}};

% Plot data color-coded according to time
uniqueTimePoint = sort(unique(timePoint));
timePointSize = zeros(numel(uniqueTimePoint), 1);
for i = 1:numel(uniqueTimePoint)
    timePointSize(i) = sum(double(timePoint == uniqueTimePoint(i)));
end


cmap = jet(numel(uniqueTimePoint));

legendLabels = cell(numel(uniqueTimePoint), 1);
legendMarkers = cell(numel(uniqueTimePoint), 1);
for i = 1:numel(legendLabels)
    minutes = uniqueTimePoint(i);
    if minutes == 0
        legendLabels{i} = sprintf('Time not specified (%d points)', timePointSize(i));
    else
        if rem(minutes, 60) == 0
            hrs = minutes / 60;
            legendLabels{i} = sprintf('%d hr (%d points)', hrs, timePointSize(i));
        else
            legendLabels{i} = sprintf('%d min (%d points)', minutes, timePointSize(i));
        end
    end
    legendMarkers{i} = 's';
end

% Create plot legend for color-coded time points
figureInitPlot(-1, -1, legendMarkers, cmap, legendLabels);
hold on;


% Plot points color coded according to time
for i = 1:numel(uniqueTimePoint)
    for j = 1:4
        markerSpec = marker{j};
        indices = (timePoint == uniqueTimePoint(i)) & (objectCount == j);
        plot(angle(indices), radius(indices), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
    end
end
title({'Polar Coordinates of Intron Locations Rotated to 1st Quadrant'; normalizationStr; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Degrees from Major Axis of Nucleus');
ylabel('Normalized Distance from Nucleus Centroid');
xlim([0, 90]);
if radiusNormalization == 1
    ylim([0, 1]);
else
    ylim([0, roundN(max(radius(:)), 1)]);
end







angleRadians = angle * (pi / 180);
X = radius .* cos(angleRadians);
Y = radius .* sin(angleRadians);

% Create plot legend for color-coded time points
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;


% Plot points color coded according to time
for i = 1:numel(uniqueTimePoint)
    for j = 1:4
        markerSpec = marker{j};
        indices = (timePoint == uniqueTimePoint(i)) & (objectCount == j);
        plot(X(indices), Y(indices), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
    end
end
title({'Cartesian Coordinates of Intron Locations Rotated to 1st Quadrant'; normalizationStr; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([0, 1]);
    ylim([0, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([0, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([0, extremeY]);
end
axis equal;





% Cluster according to polar coordinates
normalizedAngle = angle / 90;
normalizedPolar = [radius, normalizedAngle];


numPoints = numel(radius);
maxK = round(sqrt(numPoints));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints, maxK);

bestS = -Inf;
bestK = [];
bestIdx = [];
for k = 2:maxK
    [idx C] = kmeans(normalizedPolar, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(normalizedPolar, idx));
    fprintf('k: %d  mean silhouette: %f\n', k, s);
    if s > bestS
        bestS = s;
        bestK = k;
        bestIdx = idx;
    end
end

fprintf('best k between %d and %d: %d (mean silhouette: %f)\n', 2, maxK, bestK, bestS);



clusterSize = zeros(bestK, 1);
for i = 1:bestK
    clusterSize(i) = sum(double(bestIdx == i));
end

fprintf('\n');

cmap = jet(bestK);
legendLabels = cell(bestK, 1);
legendMarkers = cell(bestK, 1);
for i = 1:bestK
    legendLabels{i} = sprintf('Cluster %d (%d points)', i, clusterSize(i));
    legendMarkers{i} = 's';
end
figureInitPlot(-1, -1, legendMarkers, cmap, legendLabels);
hold on;

for i = 1:4
    markerSpec = marker{i};
    for j = 1:bestK
        indices = (objectCount == i) & (bestIdx == j);
        plot(angle(indices), radius(indices), markerSpec{1}, ...
            'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), ...
            'MarkerEdgeColor', cmap(j, :));
    end
end

title({'Polar Coordinates of Intron Locations Shifted to 1st Quadrant'; normalizationStr; 'Clustered According to Polar Coordinates'; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Degrees from Major Axis of Nucleus');
ylabel('Normalized Distance from Nucleus Centroid');
xlim([0, 90]);
if radiusNormalization == 1
    ylim([0, 1]);
else
    ylim([0, roundN(max(radius(:)), 1)]);
end






figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;
for i = 1:4
    markerSpec = marker{i};
    for j = 1:bestK
        indices = (objectCount == i) & (bestIdx == j);
        plot(X(indices), Y(indices), markerSpec{1}, ...
            'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), ...
            'MarkerEdgeColor', cmap(j, :));
    end
end

title({'Cartesian Coordinates of Intron Locations Shifted to 1st Quadrant'; normalizationStr; 'Clustered According to Polar Coordinates'; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([0, 1]);
    ylim([0, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([0, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([0, extremeY]);
end
axis equal;


numPoints = numel(radius);
maxK = round(sqrt(numPoints));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints, maxK);

cartesianPoints = [X, Y];
bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
bestDistanceTable = [];
for k = 2:maxK
    [idx C sumd D] = kmeans(cartesianPoints, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(cartesianPoints, idx));
    fprintf('k: %d   mean silhoutte: %f\n', k, s);
    if s > bestS
        bestS = s;
        bestK = k;
        bestIdx = idx;
        bestCentroids = C;
        bestDistanceTable = D;
    end
end


clusterSize = zeros(bestK, 1);
for i = 1:bestK
   clusterSize(i) = sum(double(bestIdx == i)); 
end

cmap = jet(bestK);
legendLabels = cell(bestK, 1);
legendMarkers = cell(bestK, 1);
for i = 1:bestK
    legendLabels{i} = sprintf('Cluster %d (%d points)', i, clusterSize(i));
    legendMarkers{i} = 's';
end


figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
hold on;
for i = 1:4
    markerSpec = marker{i};
    for j = 1:bestK
        indices = (objectCount == i) & (bestIdx == j);
        plot(X(indices), Y(indices), markerSpec{1}, ...
            'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), ...
            'MarkerEdgeColor', cmap(j, :));
    end
end

title({'Cartesian Coordinates of Intron Locations Shifted to 1st Quadrant'; normalizationStr; 'Clustered According to Cartesian Coordinates'; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Normalized X');
ylabel('Normalized Y');
if radiusNormalization == 1
    xlim([0, 1]);
    ylim([0, 1]);
else
    extremeX = roundN(max(abs(X(:))), 1);
    xlim([0, extremeX]);
    extremeY = roundN(max(abs(Y(:))), 1);
    ylim([0, extremeY]);
end
axis equal;


end


