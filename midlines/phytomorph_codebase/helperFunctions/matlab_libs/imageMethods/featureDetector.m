function [pointSet] = featureDetector(FileList,para,disp)     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this will be the general purpose feature detection code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           FileList    := java source of images. get will retrieve the
    %                          image from the source
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           pointSet    := List with coordinates of corners
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if disp
    if disp;wbar = waitbar(0,'Extracting features...');end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = FileList.size();
    for s = 1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load the image into qdisk
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        I = double(imread(FileList.get(s-1)));
        %%% NEEDED. SWITCH STATEMENT FOR DIFFERENT FEATURE DETECTION
        pointSet{s} = cornerFeature(I,para);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Done with feature detection:' num2str(s) ':' num2str(N) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp;wbar = waitbar(s/N,wbar,'Extracting features...');end
    end
    % if disp
    if disp;close(wbar);end
end