function cluster0(fileName)

% radiusNormalization - 1: radius divided by distance of centroid to 
% boundary; 2: radius divided by nucleus major axis length
radiusNormalization = 1;

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

% radiusColumnName = 'Marker Distance to Nucleus Centroid';    %'Normalized Radius';
% radiusColumn = find(strcmp(varnamesCA, radiusColumnName));
% assert (numel(radiusColumn) == 1, '%d occurrences of ''%s'' found', numel(radiusColumn), radiusColumnName);

angleColumnName = 'Marker Angle from Nucleus Major Axis';
angleColumn = find(strcmp(varnamesCA, angleColumnName));
assert (numel(angleColumn) == 1, '%d occurrences of ''%s'' found', numel(angleColumn), angleColumnName);

timeColumnName = 'Time in Minutes';
timeColumn = find(strcmp(varnamesCA, timeColumnName));
assert (numel(timeColumn) == 1, '%d occurrences of ''%s'' found', numel(timeColumn), timeColumnName);


[radius normalizationStr] = getMarkerDistance(varnamesCA, data, radiusNormalization);
angle = data(:, angleColumn);
objectCount = data(:, objectCountColumn);
timePoint = data(:, timeColumn);

% Shapes and sizes of plot markers
marker = {{'p', 3}, {'o', 3}, {'^', 3}, {'s', 3}};

% Plot data color-coded according to time
uniqueTimePoint = sort(unique(timePoint));
numUniqueTimePoints = numel(uniqueTimePoint);
timePointCount = zeros(numUniqueTimePoints, 1);
for i = 1:numUniqueTimePoints
   timePointCount(i) = sum(double(timePoint == uniqueTimePoint(i))); 
end

cmap = jet(numUniqueTimePoints);



legendLabels = cell(numUniqueTimePoints, 1);
legendMarkers = cell(numUniqueTimePoints, 1);
for i = 1:numel(legendLabels)
    minutes = uniqueTimePoint(i);
    if minutes == 0
        legendLabels{i} = sprintf('Time not specified  (%d points)', timePointCount(i));
    else
        if rem(minutes, 60) == 0
            hrs = minutes / 60;
            legendLabels{i} = sprintf('%d hr  (%d points)', hrs, timePointCount(i));
        else
            legendLabels{i} = sprintf('%d min  (%d points)', minutes, timePointCount(i));
        end
    end
    legendMarkers{i} = 's';
end

% Create plot legend for color-coded time points
figureInitPlot(-1, -1, legendMarkers, cmap, legendLabels);
hold on;


% Plot points color coded according to time
for i = 1:numUniqueTimePoints
   for j = 1:4
       markerSpec = marker{j};
       indices = (timePoint == uniqueTimePoint(i)) & (objectCount == j);
       plot(angle(indices), radius(indices), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
   end
end
title({'Polar Coordinates of Intron Locations'; normalizationStr; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
xlabel('Degrees from Major Axis of Nucleus');
ylabel('Normalized Distance from Nucleus Centroid');
xlim([0, 360]);
if radiusNormalization == 1
    ylim([0, 1]);
else
    ylim([0, roundN(max(radius(:)), 1)]);
end


% Plot points in Cartesian coordinates

angleRadians = angle * (pi / 180);
figureInitPlot(-2, -2, legendMarkers, cmap, legendLabels);
X = radius .* cos(angleRadians);
Y = radius .* sin(angleRadians);
hold on;

for i = 1:numUniqueTimePoints
   for j = 1:4
       markerSpec = marker{j};
       indices = (timePoint == uniqueTimePoint(i)) & (objectCount == j);
       plot(X(indices), Y(indices), markerSpec{1}, 'MarkerSize', markerSpec{2}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
   end
end
title({'Cartesian Coordinates of Intron Locations'; normalizationStr; 'Transcription Sites per Nucleus - star: 1, circle: 2, triangle: 3, square: 4'});
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
