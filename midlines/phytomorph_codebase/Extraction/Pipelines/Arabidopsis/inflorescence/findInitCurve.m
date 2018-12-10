function [curveIDX] = findInitCurve(fileList,idxL,PATHS,startPoints)
    try
        fidx = find(idxL(:,2)==1);
        I = imread(fileList{1});
        MSK = zeros(size(I));
        MSK(idxL(fidx,1)) = 1;
        MSK(:,startPoints(1,2):end) = 0;
        MSK = bwareaopen(MSK,100);
        R = regionprops(logical(MSK),'Image','PixelIdxList','Centroid');
        for h = 1:numel(R)
            value(h) = R(h).Centroid(2);
        end
        [~,sidx] = sort(value);
        R = R(sidx);
        thresh = 9;
        for e = 1:numel(R)
            MSK = zeros(size(I));
            MSK(R(e).PixelIdxList) = 1;
            vec = sum(MSK,1);
            svec = imfilter(vec,fspecial('average',[1 30]),'replicate');
            gvec = gradient(gradient(svec));
            sgvec = imfilter(gvec,fspecial('average',[1 60]),'replicate');
            
            minIDX = find((imerode(sgvec,ones(1,100)) == sgvec) & sgvec ~= 0 );
            
            maxIDX = find((imdilate(sgvec,ones(1,100)) == sgvec) & sgvec ~= 0 );
            
            fidx = find(maxIDX > minIDX(1));
            
            bidx(e) = maxIDX(fidx(1));
            %{
            zc = find(sgvec(1:end-1) < 0 & sgvec(2:end) > 0);
            bidx(e) = zc(1);
            %}
            %{
            [~,bidx(e)] = min(sgvec);
            
            
            bin = (vec < thresh) & vec ~= 0;
            bin = bwareaopen(bin,50);
            fidx = find(bin);
            bidx(e) = fidx(1);
            iL(e) = sum(bin);
            %}
            sig = mean(MSK,2);
            [~,r(e)] = max(sig);
            
            %{
            out = flattenMaskOverlay(I,logical(MSK));
            imshow(out,[]);
            hold on
            plot(bidx(e)*ones(size(MSK,1),1),1:size(MSK,1)')
            waitforbuttonpress
            %}
        end

        
        

        LEN = [];
        for sp = 1:numel(PATHS)
            t=1;
            delta = [];
            for pth = 1:numel(PATHS{sp}{t})
                tmp = PATHS{sp}{t}{pth};
                %{
                delta(pth) = norm(tmp(1,:) - [bidx(sp) r(sp)]);
                
                
               
                %}
                delta(pth) = abs(tmp(1,2) - bidx(sp));
                
                
            end
            [~,curveIDX(sp)] = min(delta);
            
            imshow(I,[]);
            hold on
            plot(PATHS{sp}{t}{curveIDX(sp)}(:,2),PATHS{sp}{t}{curveIDX(sp)}(:,1),'r');
            plot(tmp(:,2),tmp(:,1),'r')
            plot(tmp(1,2),tmp(1,2),'b*')
            plot(bidx(sp),r(sp),'g*');
            hold off
            waitforbuttonpress
            
        end
        
        
        
        
    catch ME
        ME
    end
    
    
end