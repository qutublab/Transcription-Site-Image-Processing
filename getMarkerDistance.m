function [dist normalizationStr] = getNormalizedMarkerDistance(varnamesCA, data, normalization)

centroidDistanceColumnName = 'Marker Distance to Nucleus Centroid';
centroidDistanceColumn = find(strcmp(varnamesCA, centroidDistanceColumnName));
assert (numel(centroidDistanceColumn) == 1, ...
    '%d occurrences of ''%s'' found', numel(centroidDistanceColumn), ...
    centroidDistanceColumnName);
        
centroidDist = data(:, centroidDistanceColumn);

switch normalization
    case 1
        % Divide distance from nucleus centroid to marker by length of line
        % segment from nucleus centroid to nucleus border that passes
        % through marker point
        boundaryDistanceColumnName = 'Marker Distance to Nucleus Boundary';
        boundaryDistanceColumn = find(strcmp(varnamesCA, boundaryDistanceColumnName));
        assert (numel(boundaryDistanceColumn) == 1, ...
            '%d occurrences of ''%s'' found', ...
            numel(boundaryDistanceColumn), boundaryDistanceColumnName);
        boundaryDist = data(:, boundaryDistanceColumn);
        dist = centroidDist ./ (centroidDist + boundaryDist);
        normalizationStr = 'Boundary Normalized';
    case 2
        % Divide distance from nucleus centroid to marker by length of
        % major axis of nucleus.
        nucleusMajorAxisLengthColumnName = 'Nucleus Major Axis Length';
        nucleusMajorAxisLengthColumn = find(strcmp(varnamesCA, nucleusMajorAxisLengthColumnName));
        assert (numel(nucleusMajorAxisLengthColumn) == 1, ...
            '%d occurrences of ''%s'' found', ...
            numel(nucleusMajorAxisLengthColumn), nucleusMajorAxisLengthColumnName);
        nucleusMajorAxisLength = data(:, nucleusMajorAxisLengthColumn);
        dist = centroidDist ./ nucleusMajorAxisLength;
        normalizationStr = 'Axis Normalized';
    otherwise
        error('[getMarkerDistance] Unexpected normalization parameter: %d', normalization);
end
end