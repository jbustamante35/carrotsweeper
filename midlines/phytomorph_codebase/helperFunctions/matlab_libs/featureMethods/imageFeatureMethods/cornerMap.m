function [cM] = cornerMap(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute corner map
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
    % filter the image
    I = imfilter(I,fspecial('gaussian',[2*para.vars.sig.value 2*para.vars.sig.value],para.vars.sig.value),'replicate');
    % take spactio-gradient
    D = diffmethod(I,para.vars.gradPara);                                       % get the information for the structure tensor
    % create structor tensor
    D = structureTensor(D,fspecial('gaussian',[12 12],10));                     % make the structure tensor
    % obtain the corner
    D = cornerstrength(D);                                                      % measure the tensor    
    % return corner map
    cM = D(:,:,end);
end