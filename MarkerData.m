
classdef MarkerData < handle
    properties
        fileName
        markerName
        nucleusIndex
        markerObjectIndex
        totalMarkerCount
        area
        centroidMarker
        centroidNucleus
%         normalizedRTheta
        distanceToCentroid
        distanceToBorder
        angleFromMajorAxis
        meanIntensity
        timePoint
        nucleusMajorAxisLength
        nucleusMajorAxisAngle
    end
    
    methods
        function md = MarkerData(fileName, markerName, nucleusIndex, markerObjectIndex)
            md.fileName = fileName;
            md.markerName = markerName;
            md.nucleusIndex = nucleusIndex;
            md.markerObjectIndex = markerObjectIndex;
        end
        
        function printCSV(md)
            
        end
    end
    
end
