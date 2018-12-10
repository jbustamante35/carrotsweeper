function [sI BOX] = myCrop(I,scale,BOX)
    %%%%%%%%%%
    %
    %%%%%%%%%%
    if nargin >= 2
        if scale ~= 1        
            I = imresize(I,scale);
        end
    end
    
    
    %%%%%%%%%%
    %
    %%%%%%%%%%
    if nargin >= 3
        [sI BOX] = imcrop(I,BOX);
    else
        [sI BOX] = imcrop(I);
    end
    
    
    
    %%%%%%%%%%
    %
    %%%%%%%%%%
    if nargin >= 2
        BOX = BOX*scale^-1;
    end
    
    close all;
end