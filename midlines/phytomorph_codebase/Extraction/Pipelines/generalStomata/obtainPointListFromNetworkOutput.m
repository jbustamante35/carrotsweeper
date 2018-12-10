function [pointList] = obtainPointListFromNetworkOutput(netOutput)
    % obtain the vectors
    pointList = reshape(netOutput(1,1:8),[2 4])';
end