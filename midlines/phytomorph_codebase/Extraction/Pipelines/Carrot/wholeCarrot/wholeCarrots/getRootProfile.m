function [profile,profileN] = getRootProfile(root,image_name)
    profile = [];
    profileN = [];
    % check to see if there is a carrot
    R = regionprops(any(root,2),'Area','PixelIdxList');
    % if there are objects -> carrot
    if numel(R) ~= 0 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the object with the max area and fill in only the root
        % profile
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [J,midx] = max([R.Area]);
        mtmp = zeros(size(root,1),1);
        mtmp(R(midx).PixelIdxList) = 1;
        mtmp = mtmp.*sum(root,2);
        profile = mtmp;
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trace out plongest path
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = regionprops(logical(root),'Area','PixelIdxList');
        [J,midx] = max([R.Area]);
        rootTmp = zeros(size(root));
        rootTmp(R(midx).PixelIdxList) = 1;
        toPad = 2*sum(rootTmp(1,:));
        rootTmp = [repmat(rootTmp(1,:),[toPad 1]);rootTmp];
        rootTmp = imfill(rootTmp,'holes');
        skel = bwmorph(logical(rootTmp),'skeleton',inf);
        skel(1:toPad,:) = [];
        bidx = find(imdilate(bwmorph(skel,'branchpoints',1),ones(3,3)));
        
        %skel(bidx) = 0;
        skel = bwareaopen(skel,100);
        skel = imclose(skel,ones(5));
        skel = bwlarge(skel);
        %skel = bwmorph(skel,'skeleton',inf);
        
        dist = bwdist(~logical(root));
        [x y] = find(skel);
        [idxS] = find(skel);
        
        median(dist(idxS))
        
        
        
        % construct adjacency
        DP = [x y]';
        T = Radjacency(DP,3);
        
        
        
        ep = find(bwmorph(skel,'endpoints',1));
        [ep1 ep2] = ind2sub(size(skel),ep);
        
        
        
        startIdx = find(ep1==1);
        endIdx = find(ep1~=1);
        % single trace    
        [idx(1)] = snapTo(DP',[ep1(startIdx) ep2(startIdx)]);
        [idx(2)] = snapTo(DP',[ep1(endIdx) ep2(endIdx)]);
        [path , pathcost]  = dijkstra(T , idx(2) , idx(1));
        path = DP(:,path);
        pidx = sub2ind(size(skel),path(1,:),path(2,:));
        profileN = dist(pidx);
        %}
        profileN = profile;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save root image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save the root image
        image(255*uint8(root));
        colormap(gray)
        % trace boundary
        dB = bwboundaries(logical(root));
        % display 
        hold on
        plot(dB{1}(:,2),dB{1}(:,1),'g');
        %plot(path(2,:),path(1,:),'r');
        axis off
        saveas(gca,image_name);
        close all;
    end
end