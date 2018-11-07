function [cornerMap] = cornerFeature(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % detect all corners
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script
    %                   := para{1}.value -> disk radius for nonmaximal suppresion
    %                   := para{2}.value -> corner threshold 
    %                   := para{3}.value -> sigma for gaussian filter
    %                   := para{5}.value -> mask value for clearing edges
    %                   := para{6}.value -> display during operations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           P       := List of coordinates of corners
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set the variables
    rad =       para{1}.value;          % = radius of disks
    thres =     para{2}.value;          % = threshold for minimum brightness of a corner
    sig =       para{3}.value;          % = sigma of "gaussian filtering"
    gradPara =  para{4};                % = gradient var
    mskValue =  para{5}.value;          % = mask value
    disp =      para{6}.value;          % = display (1 = on, 0 = off)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if display
    if disp;dh = figure;Iraw = I;end
    % filter the image
    I = imfilter(I,fspecial('gaussian',[2*sig 2*sig],sig),'replicate');
    % take spactio-gradient
    D = diffmethod(I,gradPara);                                 % get the information for the structure tensor
    % create structor tensor
    D = structureTensor(D,fspecial('gaussian',[12 12],10));     % make the structure tensor
    % obtain the corner
    D = cornerstrength(D,fspecial('average',[1 1]));            % measure the tensor    
    %
    cornerMap = D(:,:,end);

    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if display
    if disp;
            figure(dh);imshow(Iraw,[]);
            hold on;
            plot(P(:,2),P(:,1),'r.');
            waitforbuttonpress;
            close(dh);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
end