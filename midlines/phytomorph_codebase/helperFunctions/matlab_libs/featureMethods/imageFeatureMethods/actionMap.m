function [aM] = actionMap(fileList,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute action map over time sequence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script         
    %                   := para.sig         -> sigma for gaussian filter
    %                   := para.gradPara    -> gradient configure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           cM      := corner map       -> corner strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % preallocate memory
    [stack sz] = allocate(fileList,para);
    
    % loop over files
    N = fileList.size();
    for i = 1:N
        I = imread(fileList.get(i-1));
        if para.resize.compute
            stack(:,:,i) = imresize(I,para.resize.value);
        end
    end
    
    % calculate action map
    aM = std(stack,1,3);
    if para.resize.compute
        aM  = imresize(aM,sz);
    end
    
    % filter the action map
    h = fspecial('disk',9);
    aM = imfilter(aM,h);
    
    % normalize the action map
    aM = bindVec(aM);
end

%%%%%%%%%%%%%%%%%%
% allocate memory
function [S sz] = allocate(fileList,para)
    % peek at first image
    I = imread(fileList.get(0));
    sz = size(I);
    if para.resize.compute
        I = imresize(I,para.resize.value);
    end
    S = zeros([size(I) fileList.size()]);
end