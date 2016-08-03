function cluster10(fileName)

% radiusNormalization - 1: radius divided by distance of centroid to 
% boundary; 2: radius divided by nucleus major axis length
radiusNormalization = 1;

% Replicates parameter for kmeans function
kmeansReplicates = 1000;

locationsPerNucleus = 4;

if nargin == 0
    fileName = 'out.csv';
end

[data varnames casenames] = tblread(fileName, ',');
varnamesCA = cell(size(varnames, 1), 1);
for i = 1:numel(varnamesCA)
   varnamesCA{i} = strtrim(varnames(i, :));
end

nucleusIndexColumnName = 'Nucleus Index';
nucleusIndexColumn = find(strcmp(varnamesCA, nucleusIndexColumnName));
assert (numel(nucleusIndexColumn) == 1, '%d occurrences of ''%s'' found', numel(nucleusIndexColumn), nucleusIndexColumnName);

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
nucleusIndex = data(:, nucleusIndexColumn);

for i = 1:max(objectCount)
    count = numel(find(objectCount == i));
    nucleiCount = count / i;
    fprintf('%d nuclei with %d locations; %d total locations\n', nucleiCount, i, count);
end
fprintf('\n');

numOriginalPoints = numel(radius);



% Convert angles to radians
angle = angle * (pi / 180);

% Since two points, one at 0 degrees one at 359 degrees are near each other
% switch to cartesian coordinates to refelect this for the Matlab kmeans
% function
X = radius .* cos(angle);
Y = radius .* sin(angle);
points = [X, Y];


marker = {{'p', 3}, {'o', 3}, {'^', 3}, {'s', 3}};


% Cluster only those locations that occur in groups of a specified size within
% a nucleus
groupIndices = objectCount == locationsPerNucleus;
points = points(groupIndices, :);
numPoints = size(points, 1);
timePoint = timePoint(groupIndices);
nucleusIndex = nucleusIndex(groupIndices);
intronFileNames = casenames(groupIndices, :);



% Sanity check that every group of consecutive points comes from the
% same nucleus
for i = 1:locationsPerNucleus:numPoints
    ifn = intronFileNames(i, :);
    nucIndex = nucleusIndex(i);
    for j = (i+1):(i+locationsPerNucleus-1)
        if ~(strcmp(ifn, intronFileNames(j, :)) && nucIndex == nucleusIndex(j))
            ferror('[cluster10] Unexpected intron location sequence');
        end
    end
end

% Build table of results of comparing the points of each nucleus to the points
% of all other nuclei.
numNuclei = numPoints / locationsPerNucleus;
nucleusScores = zeros(numNuclei);
for i = 1:locationsPerNucleus:numPoints
    nucleusPoints1 = points(i:(i+locationsPerNucleus-1), :);
    tableIdx1 = ((i - 1) / locationsPerNucleus) + 1;
    for j = (i+locationsPerNucleus):locationsPerNucleus:numPoints
        nucleusPoints2 = points(j:(j+locationsPerNucleus-1), :);
        m = distanceMetric(nucleusPoints1, nucleusPoints2);
        tableIdx2 = ((j - 1) / locationsPerNucleus) + 1;
        nucleusScores(tableIdx1, tableIdx2) = m;
        nucleusScores(tableIdx2, tableIdx1) = m;
    end
end

% Cluster using kmeans
maxK = ceil(sqrt(numNuclei));
fprintf('%d data points. Begin looking for best number of clusters between 2 and %d\n', numNuclei, maxK);


bestS = -Inf;
bestK = [];
bestIdx = [];
bestCentroids = [];
bestDistanceTable = [];
for k = 2:maxK
    [idx C sumd D] = kmeans(nucleusScores, k, 'MaxIter', 10000, 'Replicates', kmeansReplicates);
    s = mean(silhouette(nucleusScores, idx));
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


% Use PCA to plot scores in 3 dimensions
[coeff, score, latent, tsquared, explained, mu] = pca(nucleusScores);
fprintf('Variance captured by first 3 principal components:');
var3 = 0;
for i = 1:3
    var3 = var3 + explained(i);
    fprintf(' %f', explained(i));
end
fprintf('   total: %f\n', var3);

pc1 = coeff(:, 1);
pc2 = coeff(:, 2);
pc3 = coeff(:, 3);


cmap = jet(bestK);

figure
for i = 1:bestK
    inCluster = bestIdx == i;
    clusterX = pc1(inCluster);
    clusterY = pc2(inCluster);
    clusterZ = pc2(inCluster);
    scatter3(clusterX, clusterY, clusterZ, 'MarkerFaceColor', cmap(i, :));
    hold on
end
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title({sprintf('Principle Components of Nuclei with %d Intron Locations', locationsPerNucleus), 'Color Coded by Cluster'});

% Plot PC2 and PC3 only
figure
for i = 1:bestK
    inCluster = bestIdx == i;
    clusterX = pc2(inCluster);
    clusterY = pc3(inCluster);
    scatter(clusterX, clusterY, 'MarkerFaceColor', cmap(i, :));
    hold on
end
xlabel('PC2');
ylabel('PC3');
title({sprintf('2nd and 3rd Principle Components of Nuclei with %d Intron Locations', locationsPerNucleus), 'Color Coded by Cluster'});


% Plot PC1 and PC3 only
figure
for i = 1:bestK
    inCluster = bestIdx == i;
    clusterX = pc1(inCluster);
    clusterY = pc3(inCluster);
    scatter(clusterX, clusterY, 'MarkerFaceColor', cmap(i, :));
    hold on
end
xlabel('PC1');
ylabel('PC3');
title({sprintf('1st and 3rd Principle Components of Nuclei with %d Intron Locations', locationsPerNucleus), 'Color Coded by Cluster'});

% Plot PC1 and PC2 only
figure
for i = 1:bestK
    inCluster = bestIdx == i;
    clusterX = pc1(inCluster);
    clusterY = pc2(inCluster);
    scatter(clusterX, clusterY, 'MarkerFaceColor', cmap(i, :));
    hold on
end
xlabel('PC1');
ylabel('PC2');
title({sprintf('1st and 2nd Principle Components of Nuclei with %d Intron Locations', locationsPerNucleus), 'Color Coded by Cluster'});



% For each cluster, plot the original intron locations
for i = 1:bestK
    inCluster = bestIdx == i;
    inClusterIndices = find(inCluster);
    clusterSize = numel(inClusterIndices);
    cmap = jet(clusterSize);
    clusterPoints = points(inCluster, :);
    figure;
    for j = 1:clusterSize
        % The values in inClusterIndices are indices in nucleusScores of
        % nuclei in the same cluster. From the inClusterIndices values,
        % recover the indices of the intron locations in each nucleus.
        firstPointIndex = ((inClusterIndices(j) - 1) * locationsPerNucleus) + 1;
        lastPointIndex = firstPointIndex + locationsPerNucleus - 1;
        % Include the first point twice so that a line will be drawn from the
        % last point to the first point
        nucleusPointsX = points([firstPointIndex:lastPointIndex, firstPointIndex], 1);
        nucleusPointsY = points([firstPointIndex:lastPointIndex, firstPointIndex], 2);
        plot(nucleusPointsX, nucleusPointsY, '-o', 'Color', [0.9 0.9 0.9], 'MarkerFaceColor', cmap(j, :), 'LineWidth', 0.01);
        hold on;
    end
    axis equal
    xlim([-1, 1]);
    ylim([-1, 1]);
    xlabel('Normalized X');
    ylabel('Normalized Y');
    title({sprintf('Cluster %d of Nuclei with %d Intron Locations', i, locationsPerNucleus), 'Color Coded by Nucleus'});

end


end
