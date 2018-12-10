function [CROPBOX] = generateCropBox(I)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for finding peaks in 1-d signals
    fS2 = 31;                                       % dilate amount for 1-d signal 
    fS1 = 201;                                      % dilate amount for 1-d signal 
    fS = 31;                                        % smoothing for measure on finding (A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for smoothing the gradient
    gradDiskSize = 11;                  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for finding the cap
    rowSampleWidth = 20;                            % thickness for sampling along rows and finding the cap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop box values
    CROPBOX_HALF_WIDTH = 100;
    CROPBOX_FORWARD = 400;
    CROPBOX_BACKWARD = 80;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain background via close operation while preserving the shape
    % of the background
    BK = imclose(double(I),strel('disk',51));
    I = I - BK;
    I = bindVec(I);        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare the edge image
    [d1 d2] = gradient(imfilter(I,fspecial('disk',gradDiskSize),'replicate'));
    G = (d1.*d1 + d2.*d2).^.5;
    G = bindVec(G);
    thresholdG = graythresh(G);
    E = G > thresholdG;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prepare for analysis
    Ii = abs(I-1);                                      % invert image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find rows for maize kernels
    s2 = sum(Ii,2).*std(abs(d1),1,2);                   % high inverted image and high gradient along x
    s2 = bindVec(s2);                                   % normalize
    s2 = imfilter(s2,fspecial('disk',fS),'replicate');  % smooth (A)
    es2 = imdilate(s2,ones(fS2,1));                     % dilate for non-max suppression
    p2 = s2 == es2;                                     % find local max
    fidx = find(p2);                                    % find local max
    sam2 = s2(fidx);                                    % sample local max
    thresh2 = graythresh(s2);                           % perform global threshold        
    thresh2 = .2;                                       % HARDWIRE THRESHOLD
    fidx = fidx(sam2 > thresh2);                        % find local max above threshold
    p2 = zeros(size(p2));                               % create zeros mask
    p2(fidx) = 1;                                       % flag local max above threshold
    P2 = repmat(p2,[1 size(I,2)]);                      % repmat mask for checker board intersection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find cols
    s1 = sum(Ii,1).*std(abs(d2),1,1);        
    s1 = bindVec(s1);
    s1 = imfilter(s1,fspecial('disk',fS),'replicate');        
    es1 = imdilate(s1,ones(1,fS1));
    p1 = s1 == es1;
    fidx = find(p1);
    sam1 = s1(fidx);
    thresh1 = graythresh(s1);
    [J sidx] = sort(sam1);
    p1 = zeros(size(p1));
    p1(fidx(sidx(end))) = 1;
    P1 = repmat(p1,[size(I,1) 1]);        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND THE CAP--sample along each row
    clear cI CROPBOX;
    fidx2 = find(p2);                                       % find the rows to sample
    for r = 1:sum(p2)            
        rowSample = E(fidx2(r)-rowSampleWidth:fidx2(r)+rowSampleWidth,:);       % sample Edge of row of thickness rowSampleWidth
        rowSample = mean(rowSample,1);                      % take the mean
        capIdx = find(rowSample ~= 0);                      % find where there is not a zero
        cI{r} = [capIdx(1) + (gradDiskSize-1)/2 fidx2(r)];  % create coordinates for cap index
        % create crop box
        UL = cI{r} - [CROPBOX_BACKWARD CROPBOX_HALF_WIDTH]; % upper left
        BR = cI{r} + [CROPBOX_FORWARD CROPBOX_HALF_WIDTH];  % bottom right
        SZ = BR - UL;                                       % size
        CROPBOX{r} = [UL SZ];                               % cropbox
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIND THE TOP AND BOTTOM--sample at each center
    CENTERS = P1.*P2;
    [cy cx] = find(CENTERS);
    TOP = {};
    BOTTOM = {};
    for r = 1:numel(cx)
        colSampleUp = E(1:cy(r),cx(r)-10:cx(r)-10);
        colSampleUp = mean(colSampleUp,2);
        colSampleDown = E(cy(r):end,cx(r)-10:cx(r)-10);
        colSampleDown = mean(colSampleDown,2);
        topIDX = find(colSampleUp);
        bottomIDX = find(colSampleDown);
        TOP{r} = [cx(r) topIDX(end) - (gradDiskSize-1)/2];
        BOTTOM{r} = [cx(r) cy(r)+bottomIDX(1) + (gradDiskSize-1)/2];
    end
end