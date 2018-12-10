function [paraSuggestion,dX] = obtainParaFromPointList(pointList)
    % obtain the major and minor axis
    majorMean = .5*(pointList(3,:) - pointList(1,:));
    minorMean = .5*(pointList(4,:) - pointList(2,:));
    % three ways to displace
    

    % normalize the vectors
    for e = 1:size(pointList,1)
        LEN(e) = norm(pointList(e,:));
        npointList(e,:) = pointList(e,:)/LEN(e);
    end


    % flip the vectors as needed
    npointList(3,:) = -npointList(3,:);
    npointList(2,:) = [npointList(2,2) -npointList(2,1)];
    npointList(4,:) = -[npointList(4,2) -npointList(4,1)];
    newAngle = atan2(npointList(:,2),npointList(:,1));

    % assign the suggested guess
    paraSuggestion = [mean(newAngle) mean(LEN([1 3])) mean(LEN([2 4])) 0 0];


    % guess at displacemenet
    dX(1,:) = .5*(majorMean + minorMean);
    dX(2,:) = mean(pointList,1);
    dX(3,:) = mean(npointList,1);
end