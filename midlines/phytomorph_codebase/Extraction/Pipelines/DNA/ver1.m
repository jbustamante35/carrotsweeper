I = imread('/home/nate/Downloads/6mar15_1mar15_PWD_m5_sp4_8_rev.tif');
I = double(I)*(2^16-1)^-1;
R = I(:,:,1);
rM = R > graythresh(R);
rM = bwareaopen(rM,100);
R = regionprops(rM,'Image','BoundingBox');
%%
close all
for e = 1:numel(R)
    try


        skeletonPoints = [];
        vec = [];
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop out the color image and get the first pass at red mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop out the color image
        tmpI = imcrop(I,R(e).BoundingBox);
        sz = size(tmpI);
        % get the mask of the object under inspection
        tmpREDMASK = imcrop(rM,R(e).BoundingBox);
        tmpREDMASK = bwlarge(tmpREDMASK);





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop distance transform for contour and midline interpolation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pad with zeros
        PD = padarray(tmpREDMASK,[10 10],0,'both');
        % generate distance transform
        DIST = bwdist(~PD) - bwdist(PD);
        % scale factor
        SV = 4;
        % resize the distance transform
        DIST = imresize(DIST,SV);
        % alpha for melting
        alpha = 1.5/2;
        % melt the distance transform
        for t = 1:50
            DIST = DIST + alpha*del2(DIST);
        end
        % get the contour from the melted transform
        dC = contour(double(DIST),[0 0]);
        % get the contour from the contour matrix
        dC = dC(:,2:dC(2,1))';
        % make the new mask
        newSCM = poly2mask(dC(:,1),dC(:,2),(sz(1)+20)*SV,(sz(2)+20)*SV);
        % skeleton
        skeletonS = bwmorph(newSCM,'thin',inf);
        % find the skelton
        [skeletonPoints(:,2),skeletonPoints(:,1)] = find(skeletonS);
        skeletonPoints = (skeletonPoints * SV^-1) - 10;
        % scale and translate the contour
        dC = (dC*SV^-1)-10;
        % display after contour and midline
        figure;
        imshow(tmpI,[]);
        hold on
        plot(dC(:,1),dC(:,2),'y')
        plot(skeletonPoints(:,1),skeletonPoints(:,2),'y.')
        hold off
        %waitforbuttonpress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop and bwlarge for the large object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % look at the blue channel for the color image
        %tmpB = tmpI(:,:,3);








        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % boundary analysis to find the end points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the boundaries
        %dB = bwboundaries(tmpREDMASK);
        % filter the contour over scales
        S = 5:1:25;
        dBS = [];
        for N = 1:numel(S)
            dBS(:,:,N) = imfilter(dC,ones(S(N),1)/S(N),'circular');
        end
        %N = 5;
        %dB{1} = imfilter(dB{1},ones(N,1)/N,'circular');
        dF = diff(dBS,1,3);
        dF = mean(dF,3);

        % measure the curvature
        %K = cwtK_closed_imfilter(dB{1},{[11]});
        K = cwtK_closed_imfilter(dC,{[11]});
        K = K.K;
        cK = find(imdilate(K,strel('disk',11)) == K);
        ridx = K(cK) < 40;
        cK(ridx) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % boundary analysis to find the end points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dijkstra trace the midline
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T = Radjacency(skeletonPoints',.5);
        idx = [];
        for r = 1:numel(cK)
            idx(r) = snapTo(skeletonPoints,dC(cK(r),:));
        end
        
        
        
        [path , pathcost]  = dijkstra(T , idx(1) , idx(2));
        path = skeletonPoints(path,:);
        
        
        
        
        % display after contour and midline and dijkstra
        figure;
        imshow(tmpI,[]);
        hold on
        plot(dC(:,1),dC(:,2),'y')
        plot(skeletonPoints(:,1),skeletonPoints(:,2),'y.')
        plot(path(:,1),path(:,2),'r')
        hold off
        %waitforbuttonpress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dijkstra trace the midline
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % threshold the blue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % threshold the blue channel
        bM = tmpB > graythresh(tmpB);
        out = flattenMaskOverlay(tmpI,logical(bM),1,'b');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % threshold the blue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %}



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sub-sample
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate curvilinear-domain
        WIDTH_NUMP = 61;
        PCA_RHO = 5;
        WIDTH = 10;
        Domain = genCurvilinearDomain(path,PCA_RHO,WIDTH,WIDTH_NUMP,[]);
        % extension on one side
        dX = -diff(Domain,1,1);
        SNIP = 50;
        dX = mean(dX(1:SNIP,:,:),1);
        dNOR = sum(dX.^2,3).^-.5;
        dX = bsxfun(@times,dX,dNOR);
        LP = mean(sum(diff(path,1,1).^2,2).^.5);
        EXT = 20;
        EXT = linspace(0,EXT,EXT/LP);
        addedP = numel(EXT);
        EXT = bsxfun(@times,EXT',dX);
        EXT = bsxfun(@plus,EXT,Domain(1,:,:));
        Domain = cat(1,flipdim(EXT,1),Domain);
        % extension on one side
        Domain = flipdim(Domain,1);
        dX = -diff(Domain,1,1);
        SNIP = 50;
        dX = mean(dX(1:SNIP,:,:),1);
        dNOR = sum(dX.^2,3).^-.5;
        dX = bsxfun(@times,dX,dNOR);
        LP = mean(sum(diff(path,1,1).^2,2).^.5);
        EXT = 20;
        EXT = linspace(0,EXT,EXT/LP);
        addedP = numel(EXT);
        EXT = bsxfun(@times,EXT',dX);
        EXT = bsxfun(@plus,EXT,Domain(1,:,:));
        Domain = cat(1,flipdim(EXT,1),Domain);
        Domain = flipdim(Domain,1);
        % reshape for sub-sampling
        dsz = size(Domain);
        DomainS = reshape(Domain,[dsz(1)*dsz(2) dsz(3)]);
        DomainS = bsxfun(@plus,DomainS,fliplr(R(e).BoundingBox(1:2)));
        vec = [];
        for k = 1:size(I,3)
            vec(:,k) = ba_interp2(I(:,:,k),DomainS(:,2),DomainS(:,1));
        end
        vec = reshape(vec,[dsz(1) dsz(2) 3]);
        %waitforbuttonpress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sub-sample
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remake curve co-ordinates based on red mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vecR = vec(:,:,1);
        vecRM = vecR > graythresh(vecR);
        vecRM = bwlarge(vecRM);
        fidx = find(any(vecRM,2));
        oldN = size(path,1);
        path = squeeze(Domain(fidx(1):fidx(end),(end-1)/2,:));
        addedP = round((size(Domain,1) - size(path,1))/2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remake curve co-ordinates based on red mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flip with blue on top
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sum blue hori
        sigB = sum(vec(:,:,3),2);
        sigBM = zeros(size(vec,1),size(vec,2));

        % find blue peak
        [~,idxB] = max(sigB);
        sigBM((idxB-8):(idxB+8),:) = 1;
        %hold on
        plot(sigB,1:size(vec,1),'b');
        % sum red vert
        sigR = sum(vec(:,:,1),1);
        % make red mask
        maskStrip = vec(:,:,1) > graythresh(vec(:,:,1));
        maskStrip = bwlarge(maskStrip);
        %
        sigM = sum(maskStrip,2);
        % make green mask
        ridx = find(maskStrip);
        greenChannel = vec(:,:,2);
        dotMask = greenChannel > graythresh(greenChannel(ridx));
        greenChannelFilter = maskStrip.*greenChannel.*dotMask;
        igreenChannelFilter = sum(greenChannelFilter,2);
        dotMaskPeaks = (imdilate(igreenChannelFilter,strel('disk',11)) == igreenChannelFilter).*igreenChannelFilter;
        dotMaskPeaks = dotMaskPeaks.*any(maskStrip,2);
        %dotMask = imdilate(dotMask,strel('disk',5,0));
        %dotMask = bwareaopen(dotMask,10);
        gR = regionprops(logical(dotMaskPeaks),'Centroid');
        %{
        rmIDX = [];
        for g = 1:numel(gR)
            if gR(g).Centroid(2) > (addedP) & gR(g).Centroid(2) < (addedP + size(path,1))
                rmIDX(g) = false;
            else
                rmIDX(g) = true;
            end
        end
        %}
        %gR(logical(rmIDX)) = [];
        %plot(1:size(vec,2),sigR,'r');
        %hold off
        %waitforbuttonpress
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flip with blue on top
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






        % measure length
        dP = diff(path,1,1);
        dL = sum(sum(dP.*dP,2).^.5,2);
        dL_c = cumsum(sum(dP.*dP,2).^.5,1);


        if idxB > size(vec,1)/2
            vec = flipdim(vec,1);
            idxB = -(idxB - size(vec,1)/2) + size(vec,1)/2;
            for g = 1:numel(gR)
                gR(g).Centroid(2) = -(gR(g).Centroid(2)-size(vec,1)/2) + size(vec,1)/2;
            end
        end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make model
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        midx = find(sigM);
        WID = 50;
        Model = zeros(size(vec)+[0 100 0]);
        Model(find(sigM),(end/2-5):(end/2+5),1) = 1;
        Model((idxB-3):(idxB+3),(end/2-8):(end/2+8),3) = 1;
        for g = 1:numel(gR)
            Model((round(gR(g).Centroid(2))-1):(round(gR(g).Centroid(2))+1),(end/2-5):(end/2+5),2) = 1;
            Model((round(gR(g).Centroid(2))-1):(round(gR(g).Centroid(2))+1),(end/2-5):(end/2+5),1) = 0;
        end


        Model((midx(1)-3):(midx(1)-1),(end/2-15):(end/2+15),1) = 1;
        Model((midx(1)-3):(midx(1)-1),(end/2-15):(end/2+15),2) = 1;
        Model((midx(end)+1):(midx(end)+3),(end/2-15):(end/2+15),1) = 1;
        Model((midx(end)+1):(midx(end)+3),(end/2-15):(end/2+15),2) = 1;

        close all
        imshow(Model,[]);

        %Model = padarray(Model,[0 50],0,'both');
        imshow([vec Model],[]);

        
        
        text(size(vec,2)+size(Model,2)/2-4-30,midx(1)-2,'0','Color','y');
        text(size(vec,2)+size(Model,2)/2-4-30,midx(end)+2,num2str(round(sum(dL))),'Color','y');
        
        
        text(size(vec,2)+size(Model,2)/2-4-60,idxB,[num2str(round(dL_c(idxB-addedP+1))) '--------------' num2str(dL_c(idxB-addedP+1)/sum(dL),2)],'Color','b');
        for g = 1:numel(gR)
            text(size(vec,2)+size(vec,2)-40-WID/2+50,round(gR(g).Centroid(2)),[num2str(round(dL_c(round(gR(g).Centroid(2)-addedP+1)))) '-----' num2str((dL_c((gR(g).Centroid(2)-addedP+1)))/sum(dL),2)],'Color','g');
        end
        
        
        
        drawnow
        imwrite([vec Model],['/home/nate/Downloads/forApril2/'  num2str(e) '_straight.tif'])
        %waitforbuttonpress
        
        close all
    catch ME
        break
    end
    close all
end


%%
close all
imshow(rM)
for e = 1:numel(R)
    tmpI = double(R(e).Image);
    skel{e} = bwmorph(tmpI,'skeleton',inf);
    tmpI = cat(3,tmpI,tmpI,skel{e});
    imshow(double(tmpI),[])
    drawnow
    waitforbuttonpress
    
end