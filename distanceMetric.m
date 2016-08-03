
% Arguments pointArray1 and pointArray2 are n-by-2 arrays containing the
% cartesian coordinates of n points.  This function matches points in the two
% arrays by minimizing the aggregate difference in distance (sum of squares)
% bewteen points in pointArray1 and corresponding points in pointArray2.
% The least possible aggregate distance is returned.


function lowestSumSquares = distanceMetric(pointArray1, pointArray2)
numPoints = size(pointArray1, 1);
assert(numPoints == size(pointArray2, 1), '[distanceMetric] array arguments are of different sizes');


%fprintf('   A      B      C\n');
%for i = 1:numPoints
%    fprintf('(%d, %d) ', pointArray1(i, 1), pointArray1(i, 2));
%end
%fprintf('\n');
%for i = 1:numPoints
%    fprintf('(%d, %d) ', pointArray2(i, 1), pointArray2(i, 2));
%end
%fprintf('\n');


% Generate permutations for indices of pointArray2.  This will allow the
% examination of all possible ways of matching points in pointArray1 with
% points in pointArray2.
permutations = perms(1:numPoints);
numPermutations = size(permutations, 1);


% indexPairs holds all possible unique, nonduplicating pairs of indices,
% one pair per row.  It is used to select pairs of points for distance
% measurements.
indexPairs = nonduplicatePairs(numPoints);
numIndexPairs = size(indexPairs, 1);


lowestSumSquares = Inf;
bestPermutation = [];
for i = 1:numPermutations
    onePermutation = permutations(i, :);

%fprintf('permutation:');
%for j = 1:numel(onePermutation)
%    fprintf(' %d', onePermutation(j));
%end
%fprintf('\n');

    distanceSumSquares = 0;
    for j = 1:numIndexPairs
        p1 = indexPairs(j, 1);
        p2 = indexPairs(j, 2);
        % Compute distance between points in pointArray1
        x1 = pointArray1(p1, 1);
        y1 = pointArray1(p1, 2);
        x2 = pointArray1(p2, 1);
        y2 = pointArray1(p2, 2);
        d1 = distance(x1, y1, x2, y2);
        % Compute distance between corresponding points in pointArray2
        x1 = pointArray2(onePermutation(p1), 1);
        y1 = pointArray2(onePermutation(p1), 2);
        x2 = pointArray2(onePermutation(p2), 1);
        y2 = pointArray2(onePermutation(p2), 2);
        d2 = distance(x1, y1, x2, y2);

%fprintf('d1(%s, %s) = %f  d2(%s, %s) = %f  (d1-d2)^2 = %f\n', 'A'+(p1-1), 'A'+(p2-1), d1, 'A'+(onePermutation(p1)-1), 'A'+(onePermutation(p2)-1), d2, (d1 - d2)^2);  
        distanceSumSquares = distanceSumSquares + ((d1 - d2)^2);
    end

%fprintf('total = %f\n', distanceSumSquares);

    if distanceSumSquares < lowestSumSquares
        bestPermutation = onePermutation;
    end
    lowestSumSquares = min(lowestSumSquares, distanceSumSquares);
end
%bestPermutation    
end



% Returns an m-by-2 array containing the m possible unique pairs of n elements
% where no element is paired with itself

function pairArray = nonduplicatePairs(n)
% The number of possible unique pairs is the sum of numbers from 1 to n-1.
numPairs = (n - 1) * n / 2;
pairArray = zeros(numPairs, 2);
row = 1;
for i = 1:(n-1)
    for j = (i+1):n
        pairArray(row, :) = [i, j];
        row = row + 1;
    end
end
end

function d = distance(x1, y1, x2, y2)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end
