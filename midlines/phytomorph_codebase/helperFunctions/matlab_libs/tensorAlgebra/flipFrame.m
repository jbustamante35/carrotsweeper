function [nframe] = flipFrame(frame)   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flipFrame 2D and via (-)op on primary vector and twist 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %       frame = frame ?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       nframe = new frame ?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nframe = [-frame(:,1) twistVec(-frame(:,1))];
    if (size(frame,2) == 3)
        nframe = [nframe frame(:,3)];
    end
end
