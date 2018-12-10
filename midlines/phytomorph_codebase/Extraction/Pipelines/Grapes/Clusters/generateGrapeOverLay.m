function []  = generateGrapeOverLay(I,M,fileName)
    % clear some border
    vec = M(1,:);
    M(1,:) = 0;
    M = imclearborder(M);
    M(1,:)= vec;
    % area open
    M = bwareaopen(M,1000);
    if sum(M(:)) > 3000
        % normalize color image to 1
        for k = 1:size(I,3)
            tmp = I(:,:,k);
            if max(tmp) > 1
                tmp = tmp / 255;
            end
            I(:,:,k) = tmp;
        end
        
        
        
        % create overlay
        out = flattenMaskOverlay(I,logical(M),.5,'b');
        image(out);
        axis off
        hold on
        
        % sweep mass lines
        for e = 1:30:size(M,1)
            strip = M(e,:);
            R = regionprops(logical(strip),'PixelIdxList');
            X = find(strip);
            for r = 1:numel(R)
                plot([R(r).PixelIdxList(1) R(r).PixelIdxList(end)],[e e],'g')
            end
            %{
            if ~isempty(X)
                plot([X(1) X(end)],[e e],'g')
            end
            %}
        end
        R = regionprops(logical(M),'BoundingBox');
        for e = 1:numel(R)
            rectangle('Position',R(e).BoundingBox,'EdgeColor','r');
        end
        % get boundaries
        B = bwboundaries(M,8,'holes');
        % plot boundaries
        for b = 1:numel(B)
            plot(B{b}(:,2),B{b}(:,1),'g')
        end
        % fill holes
        hM = imfill(M,'holes');
        h = hM == 1 & M == 0;
        R = regionprops(logical(h),'Area','Centroid','EquivDiameter');
        % get centroids of holes
        centroids = [];
        for e = 1:numel(R)
            centroids = [centroids ; R(e).Centroid];
            RAD(e) = R(e).EquivDiameter;
        end
        % show holes
        if ~isempty(centroids)
            viscircles(centroids,RAD,'LineWidth',1);
        end
        drawnow
        axis equal
        % save and close
        saveas(gca,fileName);
        close all
    end
    
end