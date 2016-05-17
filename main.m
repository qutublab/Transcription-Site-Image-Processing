

function main()

% dirNameCA = {'/home/bl6/Mancini/Data/Figure 1 _ 24h', ...
%     '/home/bl6/Mancini/Data/Figure 2 _ 30min', ...
%     '/home/bl6/Mancini/Data/Figure 3 _ 30min'};

dirNameCA = {'/home/bl6/Mancini/Data/MCF-7 time course'};

dapiId = 'w435';
greb1ExonId = 'w594';
greb1IntronId = 'w676';

stats = cell(numel(dirNameCA), 1);
for i =  1:numel(dirNameCA)
    dirName = dirNameCA{i};
    %     d = dir([dirName, '/*', dapiId, '-max.tif']);
    d = dir([dirName, '/*', dapiId, '.tif']);
    statsDir = cell(numel(d), 1);
    for j = 1:numel(d)
        if strcmp(d(j).name, '20151026_2hrs_Greb_EI_5_R3D_D3D_PRJ_w435.tif')
            fprintf('[main] Skipping 20151026_2hrs_Greb_EI_5_R3D_D3D_PRJ_w435.tif due to intron segmentation problems\n');
            continue;
        end
        dapiFileName = [dirName, '/', d(j).name];
        intronFileName = replaceLast(dapiFileName, dapiId, greb1IntronId);
        timeInMinutes = extractTimeInMinutes(dapiFileName);
        D = segmentDapi(dapiFileName);
        I = segmentGreb1Intron(intronFileName);
        intronImage = mat2gray(imread(intronFileName));
        [dapiLabel numDapiLabels] = bwlabel(D);
        dapiProps = regionprops(dapiLabel, 'Centroid', 'Orientation');
        statsFile = cell(numDapiLabels, 1);
        for k = 1:numDapiLabels
            dapiMask = (dapiLabel == k);
            dapiCentroid = dapiProps(k).Centroid;
            majAxisLen = majorAxisLength(dapiMask, dapiCentroid, dapiProps(k).Orientation);
            intronMask = I & dapiMask;
            cc = bwconncomp(intronMask);
            if (cc.NumObjects > 4)
                fprintf('[main] Skipping %d Intron objects in nucleus %d in file %d %s\n', cc.NumObjects, k, j, intronFileName);
                continue;
            end
            intronCentroids = computeCentroids(cc.PixelIdxList, size(intronMask));
%             fprintf('Invoking polar.m for nucleus %d of nucleus file %d %s\n', k, i, dapiFileName);
            rThetaArr = polar(intronCentroids, dapiMask, dapiCentroid, dapiProps(k).Orientation);
            statsNucleus = cell(cc.NumObjects, 1);
            for m = 1:cc.NumObjects
                md = MarkerData(intronFileName, 'Intron', k, m);
                md.totalMarkerCount = cc.NumObjects;
                md.area = numel(cc.PixelIdxList{m});
                md.centroidMarker = intronCentroids{m};
                md.centroidNucleus = dapiCentroid;
%                 md.normalizedRTheta = rThetaArr(m, :);
                md.distanceToCentroid = rThetaArr(m, 1);
                md.distanceToBorder = rThetaArr(m, 2) - rThetaArr(m, 1); % distance from marker to nucleus border
                md.angleFromMajorAxis = rThetaArr(m, 3);
                intronPixels = intronImage(cc.PixelIdxList{m});
                md.meanIntensity = mean(intronPixels(:));
                md.timePoint = timeInMinutes;
                md.nucleusMajorAxisLength = majAxisLen;
                md.nucleusMajorAxisAngle = -dapiProps(k).Orientation;
                statsNucleus{m} = md;
            end
            statsFile{k} = statsNucleus;
        end
        statsDir{j} = statsFile;
    end
    stats{i} = statsDir;
end

outFileName = 'out.csv';
writeFileAll(stats, outFileName);
fprintf('Wrote file %s\n', outFileName);

daveOutFileName = 'outDave.csv';
writeFile4Dave(stats, daveOutFileName);
fprintf('Wrote file %s\n', daveOutFileName);

end




function s2 = replaceLast(s, old, new)
k = strfind(s, old);
if numel(k) == 0
    error('[intronPolar.replaceLast] Unable to find substring %s string %s', s, old);
end

prefix = s(1:(k(end)-1));
suffix = s((k(end)+4):end);

s2 = [prefix, new, suffix];

end



function centroidCA = computeCentroids(linearIdxList, sz)
numCentroids = numel(linearIdxList);
centroidCA = cell(numCentroids, 1);
for i = 1:numCentroids
    [R C] = ind2sub(sz, linearIdxList{i});
    centroidCA{i} = [sum(C(:)) / numel(C), sum(R(:)) / numel(R)];
end
end

function markerName = firstMarkerName(stats)
markerName = '';
while isempty(markerName)
    for i = 1:numel(stats)
        statsDir = stats{i};
        for j = 1:numel(statsDir)
            statsFile = statsDir{j};
            for k = 1:numel(statsFile)
                statsNucleus = statsFile{k};
                for m = 1:numel(statsNucleus)
                    md = statsNucleus{m};
                    markerName = md.markerName;
                    return;
                end
            end
        end
    end
end
end

function writeFile(stats, outFileName);
[fid errmsg] = fopen(outFileName, 'w');
assert(fid ~= -1, errmsg);
markerName = firstMarkerName(stats);
header = {[markerName ' File Name'], 'Time in Minutes', 'Nucleus Index', 'Total Marker Object Count', 'Marker Object Index', 'Normalized Radius', 'Angle from Nucleus Major Axis'};
for i = 1:numel(header)
    if i > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', header{i});
end
fprintf(fid, '\n');
for i = 1:numel(stats)
    statsDir = stats{i};
    for j = 1:numel(statsDir)
        statsFile = statsDir{j};
        for k = 1:numel(statsFile)
            statsNucleus = statsFile{k};
            if isempty(statsNucleus)
                continue;
            end
            for m = 1:numel(statsNucleus)
                md = statsNucleus{m};
                fprintf(fid, '%s,%d,%d,%d,%d,%f,%f\n', ...
                    md.fileName, md.timePoint, md.nucleusIndex, ...
                    md.totalMarkerCount, md.markerObjectIndex, ...
                    md.normalizedRTheta(1), md.normalizedRTheta(2));
            end
        end
    end
end
end

function writeFileAll(stats, outFileName);
[fid errmsg] = fopen(outFileName, 'w');
assert(fid ~= -1, errmsg);
markerName = firstMarkerName(stats);
header = {[markerName ' File Name'], 'Time in Minutes', ...
    'Nucleus Index', 'Nucleus Centroid X', 'Nucleus Centroid Y', ...
    'Nucleus Major Axis Length', 'Nucleus Major Axis Angle', ...
    'Total Marker Object Count', 'Marker Object Index', ...
    'Marker Centroid X', 'Marker Centroid Y', ...
    'Marker Distance to Nucleus Centroid', ...
    'Marker Distance to Nucleus Boundary', ...
    'Marker Angle from Nucleus Major Axis'};
for i = 1:numel(header)
    if i > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', header{i});
end
fprintf(fid, '\n');
for i = 1:numel(stats)
    statsDir = stats{i};
    for j = 1:numel(statsDir)
        statsFile = statsDir{j};
        for k = 1:numel(statsFile)
            statsNucleus = statsFile{k};
            if isempty(statsNucleus)
                continue;
            end
            for m = 1:numel(statsNucleus)
                md = statsNucleus{m};
                fprintf(fid, '%s,%d,%d,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f\n', ...
                    md.fileName, md.timePoint, md.nucleusIndex, ...
                    md.centroidNucleus(1), md.centroidNucleus(2), ...
                    md.nucleusMajorAxisLength, ...
                    md.nucleusMajorAxisAngle, md.totalMarkerCount, ...
                    md.markerObjectIndex, md.centroidMarker(1), ...
                    md.centroidMarker(2), md.distanceToCentroid, ...
                    md.distanceToBorder, md.angleFromMajorAxis);
            end
        end
    end
end
end

function writeFile4Dave(stats, fileName)
[fid errmsg] = fopen(fileName, 'w');
assert(fid ~= -1, errmsg);
markerName = firstMarkerName(stats);
header = cell(4 + 16, 1);
header{1} = [markerName, ' File Name'];
header{2} = 'Nucleus Index';
header{3} = 'Nucleus Centroid X';
header{4} = 'Nucleus Centroid Y';
for i = 1:4
    header{(i * 4) + 1} = sprintf('X Spot %d', i);
    header{(i * 4) + 2} = sprintf('Y Spot %d', i);
    header{(i * 4) + 3} = sprintf('Area Spot %d', i);
    header{(i * 4) + 4} = sprintf('Average Intensity Spot %d', i);
end

for i = 1:numel(header)
    if i > 1
        fprintf(fid, ',');
    end
    fprintf(fid, header{i});
end
fprintf(fid, '\n');
for i = 1:numel(stats)
    statsDir = stats{i};
    for j = 1:numel(statsDir)
        statsFile = statsDir{j};
        for k = 1:numel(statsFile)
            statsNucleus = statsFile{k};
            if numel(statsNucleus) == 0 continue; end
            md = statsNucleus{1};
            fprintf(fid, '%s,%d,%f,%f', md.fileName, md.nucleusIndex, md.centroidNucleus(1), md.centroidNucleus(2));
            for m = 1:numel(statsNucleus)
                md = statsNucleus{m};
                fprintf(fid, ',%f,%f,%d,%f', md.centroidMarker(1), md.centroidMarker(2), md.area, md.meanIntensity);
            end
            fprintf(fid, '\n');
        end
    end
end

end
