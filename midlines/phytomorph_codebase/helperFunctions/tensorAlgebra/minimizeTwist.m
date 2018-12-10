function [NewsourceFrame] = minimizeTwist(sourceFrame,targetFrame)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % minimize twist energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %       sourceFrame = source frame 
    %       targetFrame = target frame 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       NewsourceFrame = new source frame ?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    potentialFrameSet = cat(3,sourceFrame,flipFrame(sourceFrame));
    objFunc = targetFrame(:,1)'*squeeze(potentialFrameSet(:,1,:));
    [mn argmn] = max(objFunc);
    NewsourceFrame = potentialFrameSet(:,:,argmn);
end
