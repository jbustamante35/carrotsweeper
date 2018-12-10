function [MASK gBNK] = createKernelEdgeMask(fileName)

        gBNK = [];
        warning off
        GLOBAL_UPPERLEFT = [10 300];

        disp = 0;
        I = imread(fileName);        
        MASK = zeros(size(I));
        IB = I;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find all corner(s) in the first image
        cornersFeatures = cornerFeatureMap(5);
        cornerFeatures = cornersFeatures.computeFeatureMap(double(IB));
        allCorners = simplePointExtractor(20,.000008);
        allCorners = allCorners.extractPoints(cornerFeatures);
        
        
        I = double(I(GLOBAL_UPPERLEFT(1):end-10,GLOBAL_UPPERLEFT(2):end-500))/255;
        % for non max
        fS2 = 31;
        fS1 = 201;
        fS = 31;                % smoothing for measure on finding (A)
        gradDiskSize = 11;
        CROPBOX_HALF_WIDTH = 100;
        % for plots
        mag2 = 500;
        mag1 = 3000;
        % row samples width
        rowS = 20;              % thickness for sampling along rows and finding the cap
        curveBank = {};
    
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
        sam2 = s2(fidx);                                    % sampel local max
        thresh2 = graythresh(s2);                           % perform global threshold        
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
            rowSample = E(fidx2(r)-rowS:fidx2(r)+rowS,:);       % sample Edge of row of thickness rowS
            rowSample = mean(rowSample,1);                      % take the mean
            capIdx = find(rowSample ~= 0);                      % find where there is not a zero
            cI{r} = [capIdx(1) + (gradDiskSize-1)/2 fidx2(r)];  % create coordinates for cap index
            % create crop box
            UL = cI{r} - [100 CROPBOX_HALF_WIDTH];               % upper left
            BR = cI{r} + [500 CROPBOX_HALF_WIDTH];              % bottom right
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


        
        if disp
            imshow(I,[])
            hold on;
            plot(mag2*p2,1:size(I,1),'r');
            plot(1:size(I,2),mag1*p1,'g');
            %
            for i = 1:numel(cI)
                plot(cI{i}(1),cI{i}(2),'b*')
            end
            %
            for i = 1:numel(TOP)
                plot(TOP{i}(1),TOP{i}(2),'c*')
                plot(BOTTOM{i}(1),BOTTOM{i}(2),'y*')
            end
            %
            for i = 1:numel(TOP)
                plot(cx(i),cy(i),'m*')
            end
        end
        
        
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % analysis of each image box
        CUR = myHS_X('phytoAcurve');
        for i = 1:numel(CROPBOX)
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % crop and filter
            sI{i} = imcrop(I,CROPBOX{i});            
            sI{i} = imfilter(sI{i},fspecial('gaussian',11,3),'replicate');            
            C = contourc(sI{i},50);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % get the curve structure
            str = 1;
            c = 1;
            clear curve
            while str < size(C,2)
                ed = str + C(2,str);
                curve(c).level = C(1,str);
                curve(c).data = C(:,str+1:ed);
                curve(c).length = size(curve(c).data,2);
                c = c + 1;
                str = ed + 1;
            end
            [curve] = selectClosedCurves(curve);
            [J,sidx] = sort([curve.length]);
            % make mask from longest curve
            BW = poly2mask(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),size(sI{i},1),size(sI{i},2));
            % perform the distance transform on the mask
            sDIST = double(bwdist(~BW));
            
            
            
            
            %{
            SNIP = 30;
            BANK = [curve(sidx(end)).data curve(sidx(end)).data curve(sidx(end)).data];
            for e = 1:(size(BANK,2) - SNIP)
                BNK(:,:,e) = bsxfun(@minus,BANK(:,e:e+SNIP),BANK(:,e+SNIP/2));
            end
            BNK = BNK(:,:,size(curve(sidx(end)).data,2):2*size(curve(sidx(end)).data,2)-1);
            %}
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the highest curvature for the longest curve
            % first smooth large and find the biggest peak in the curvature
            % profile then find the curve within a window which is
            % contributing most the the bend
            LOC = bsxfun(@plus,curve(sidx(end)).data,CROPBOX{i}(1:2)');
            LOC = bsxfun(@plus,LOC,flipud(GLOBAL_UPPERLEFT'));
            LOC = round(LOC);
            for e = 1:size(LOC,2)
                MASK(LOC(2,e),LOC(1,e)) = 1;
            end
        end
end