function [] = main_condor(fileList,smoothValue,threshSIG,EXT,topTRIM,SNIP)
    % topTRIM = .45*size(I,1))
    % read the image
    I = imread(fileList);
    % make gray scale
    G = rgb2gray(I);
    % filter the image
    G = imfilter(G,fspecial('gaussian',[13 13],4));
    % find edge
    E = edge(G);
    % integrate
    sig = sum(E,1);
    % smooth
    sig = imfilter(sig,fspecial('average',[1 smoothValue]),'replicate');
    % find thhe gaps
    BLOCK = sig < threshSIG;
    % remove the non-gaps that are less than 50
    BLOCK = bwareaopen(BLOCK,50);
    % close the pot holder chunks
    BLOCK = imclose(BLOCK,strel('disk',100));
    % extend the conetainers holder blocks
    eBLOCK = imerode(BLOCK,strel('disk',[EXT]));
    % make an image mask
    MASK = repmat(eBLOCK,[size(I,1) 1]);
    % get the bounding boxs for each mask
    R = regionprops(~MASK,'BoundingBox');
    % for each cone-tainer
    for e = 1:numel(R)
        % crop the strip
        tmpD = imcrop(I,R(e).BoundingBox);
        % trim the top
        tmpD(round(1:topTRIM:,:) = [];
        % get the size
        SZ = size(tmpD);
        % make plant mask
        MASK = getMASK_ver0(tmpD);
        % if the mask is blank
        if sum(MASK(:))/prod(size(MASK)) < thresP
            % pad the array
            tmpMASK = padarray(MASK, [300 0], 'replicate', 'post');
            % get the skeleton
            SKEL = bwmorph(tmpMASK,'thin',inf);
            % 
            SKEL = SKEL(1:size(SKEL,1),:);
            % find the skelton for tracing
            [r c] = find(SKEL);
            % sum the mask for the height calculation
            sig = sum(MASK,2);
            % find the pixlels
            fidx = find(sig);
            % find the top pixel
            HEIGHT = fidx(1);
            % find the biomass
            dBIOMASS = sum(MASK(:));
            % find the tips
            EP = imfilter(double(SKEL),ones(3,3));
            [re ce] = find(EP == 2 & SKEL);
            % SNIP off some stem
            baseMASK = sum(MASK(end-SNIP:end,:),1);
            % mean along the stem snip
            basePoint(1) = mean(find(baseMASK));
            % set the basepoint 1 to the size of the mask
            basePoint(2) = size(MASK,1);
            % find skeleton
            [x y] = find(SKEL);
            % stack the skeleton points for tracing
            DP = [x y]';
            % make adjaceny matrix
            T = Radjacency(DP,3);
            % find the longest path from the stem end point to the leaf tip
            pathcost = [];
            path = {};
            % snap the base point to the skeleton
            [idx(1)] = snapTo(DP',[fliplr(basePoint)]);
            for i = 1:numel(re)
                % find the end point in the skeleton
                [idx(2)] = snapTo(DP',[re(i) ce(i)]);
                % trace
                [path{i} , pathcost(i)]  = dijkstra(T , idx(1) , idx(2));
            end
            % set pathcost of inf to zero
            pathcost(isinf(pathcost)) = 0;
            % find the zeros - including inf path cost
            ridx = find(pathcost==0);
            % remove the 0 length paths
            pathcost(ridx) = [];
            % remove the 0 length paths
            path(ridx) = [];
            % find the max path cost
            [J,midx] = max(pathcost);
            % make mask overlay
            out = flattenMaskOverlay(double(cSTORE{e})/255, MASK,.65,'g');
            
            fprintf(['starting with image and results display \n']);
            % display the results            
            h = image(I);
            imshow(out,[])
            hold on
            plot(1:SZ(2),HEIGHT*ones(1,SZ(2)),'g')
            %plot(c,r,'r.')
            plot(ce,re,'b*')
            for i = 1:numel(pathcost)
                plot(DP(2,path{i}),DP(1,path{i}),'r','LineWidth',3)  
            end
            plot(DP(2,path{midx}),DP(1,path{midx}),'k','LineWidth',3)
            hold off
            drawnow 
            
    end
end