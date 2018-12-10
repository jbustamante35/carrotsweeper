
function [P] = featureSelection(logicalStack)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % feature selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           logicalStack = stack of binary ("logical") images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           P = List of coordinates of corners
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % AND
    logicalResult = all(logicalStack,3);
    % find points
    [r,c] = find(logicalResult);
    % return points
    P = [c r];
end
