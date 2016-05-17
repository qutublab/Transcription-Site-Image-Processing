function cluster1(fileName)

% radiusNormalization - 1: radius divided by distance of centroid to 
% boundary; 2: radius divided by nucleus major axis length
radiusNormalization = 2;

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

% Shapes and sizes of plot markers
marker = {{'p', 3}, {'o', 3}, {'^', 3}, {'s', 3}};


% Begin cluster analysis
numPoints = numel(radius);
maxK = round(sqrt(numPoints));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numPoints, maxK);
bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
for k = 2:maxK
    [idx C] = kmeans(radius, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(radius, idx));
    fprintf('k: %d   mean silhoutte: %f\n', k, s);
    if s > bestS
        bestS = s;
        bestK = k;
        bestIdx = idx;
        bestCentroids = C;
    end
end

fprintf('Best k between %d and %d: %d (mean silhouette: %f)\n', 2, maxK, bestK, bestS);

clusterSize = zeros(bestK, 1);
for i = 1:bestK
   clusterSize(i) = sum(double(bestIdx == i)); 
end

% Plot points according to angle and radius, grouped by number of marker
% objects in nucleus.
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
    
%    indices = bestIdx == i;
    plot(angle(indices), radius(indices), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(j, :), 'MarkerEdgeColor', cmap(j, :));
    end
end

[~, sortedIndices] = sort(bestCentroids);
cluster = radius(bestIdx == sortedIndices(1));
clusterMax = max(cluster);
for i = 2:(bestK)
    cluster2 = radius(bestIdx == sortedIndices(i));
    cluster2Min = min(cluster2);
    mid = (clusterMax + cluster2Min) / 2;
    plot([0, 360], [mid, mid], '-k');
    clusterMax = max(cluster2);
end

title({'Polar Coordinates of Intron Locations'; normalizationStr; 'Clustered According to Radius'; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Degrees from Major Axis of Nucleus');
ylabel('Normalized Distance from Nucleus Centroid');
xlim([0, 360]);
if radiusNormalization == 1
    ylim([0, 1]);
else
    ylim([0, roundN(max(radius(:)), 1)]);
end


end
