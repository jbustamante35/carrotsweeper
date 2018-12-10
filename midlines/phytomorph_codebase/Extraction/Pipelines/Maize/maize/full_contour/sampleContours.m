function [data] = sampleContours(fileStack,disp,outPath)
    try
        warning off
        mkdir(outPath);
        
        GLOBAL_UPPERLEFT = [10 300];
        fileName = fileStack{1};
        I = imread(fileName);
        IB = I;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find all corner
        cornersFeatures = cornerFeatureMap(5);
        cornerFeatures = cornersFeatures.computeFeatureMap(double(IB));
        allCorners = simplePointExtractor(20,.000008);
        allCorners = allCorners.extractPoints(cornerFeatures);
        
        
        I = double(I(GLOBAL_UPPERLEFT(1):end-10,GLOBAL_UPPERLEFT(2):end-500))/255;
        % for non max
        fS2 = 31;
        fS1 = 201;
        fS = 31;
        gradDiskSize = 11;
        % for plots
        mag2 = 500;
        mag1 = 3000;
        % row samples width
        rowS = 10;
        curveBank = {};
        % for sampling curves
        T = 30;
        R = 50;
        nhSZ = [T R];
        [THETA RAD] = ndgrid(linspace(-pi,pi,T),linspace(0,30,R));
        NH = [RAD(:).*cos(THETA(:)) RAD(:).*sin(THETA(:))]';
        %{
        [S1 S2] = ndgrid(linspace(-100,100,200),linspace(-100,100,200));
        NH = [S1(:) S2(:)]';
        nhSZ = [200 200];
        %}
    
        BK = imclose(double(I),strel('disk',51));
        I = I - BK;
        I = bindVec(I);


        [d1 d2] = gradient(imfilter(I,fspecial('disk',gradDiskSize),'replicate'));
        G = (d1.*d1 + d2.*d2).^.5;
        G = bindVec(G);
        thresholdG = graythresh(G);
        E = G > thresholdG;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prepare for rows
        Ii = abs(I-1);      % invert image
        s2 = sum(Ii,2).*std(abs(d1),1,2);   % high inverted image and high gradient along x
        s2 = bindVec(s2);   % normalize
        s2 = imfilter(s2,fspecial('disk',fS),'replicate');  % smooth
        es2 = imdilate(s2,ones(fS2,1));     % dilate for non-max suppression
        p2 = s2 == es2;     % find local max
        fidx = find(p2);    % find local max
        sam2 = s2(fidx);    % sampel local max
        thresh2 = graythresh(s2);           % perform global threshold        
        fidx = fidx(sam2 > thresh2);        % find local max above threshold
        p2 = zeros(size(p2));               % create zeros mask
        p2(fidx) = 1;                       % create mask
        P2 = repmat(p2,[1 size(I,2)]);      % repmat mask

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


        %%% find plug sites
        %fM = levelSetKurvatureFeatureMap();
        %fM = fM.computeFeatureMap(imfilter(I,fspecial('disk',11),'replicate'));

        CENTERS = P1.*P2;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND THE CAP--sample along each row
        clear cI CROPBOX;
        fidx2 = find(p2);                   % find the rows to sample
        for r = 1:sum(p2)            
            rowSample = E(fidx2(r)-rowS:fidx2(r)+rowS,:);   % sample Edge of row of thickness rowS
            rowSample = mean(rowSample,1);  % take the mean  
            capIdx = find(rowSample ~= 0);  % find where there is not a zero
            cI{r} = [capIdx(1) + (gradDiskSize-1)/2 fidx2(r)]; % create coordinates for cap index
            % create crop box
            UL = cI{r} - [20 100];  % upper left
            BR = cI{r} + [500 100]; % bottom right
            SZ = BR - UL;           % size
            CROPBOX{r} = [UL SZ];   % cropbox
        end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIND THE TOP AND BOTTOM--sample at each center
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
        
        
        
        
        % analysis of each image box
        CUR = myHS_X('phytoAcurve');
        for i = 1:numel(CROPBOX)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % crop and filter
            sI{i} = imcrop(I,CROPBOX{i});            
            sI{i} = imfilter(sI{i},fspecial('gaussian',11,5),'replicate');
            hold on;
            C = contourc(sI{i},10);
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
            sDIST = double(bwdist(~BW));
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the highest curvature
            KSNIP = 50;
            out = cwtK_closed(curve(sidx(end)).data',{60});
            [J tipIDX] = min(out.K);
            out = cwtK_closed(curve(sidx(end)).data',{10});
            [J fine_tipIDX] = min(out.K(tipIDX-KSNIP:tipIDX+KSNIP));
            tipIDX = tipIDX + (fine_tipIDX - KSNIP);
            tmpPoint = (curve(sidx(end)).data(:,tipIDX)' + CROPBOX{i}(1:2) + fliplr(GLOBAL_UPPERLEFT));
            % find the nearest corner
            dist = [];
            for e = 1:numel(allCorners)
                tpt = allCorners{e}(1);
                dist(e) = norm(tmpPoint - tpt(1:2));
            end
            [J selidx] = min(dist);            
            CUR{i} = phytoAcurve.contructFromApoint(allCorners{selidx});
            
            
            
            
            
            % constuct init information for tracing midline
            d1X1 = cwt(curve(sidx(end)).data(1,:),7,'gaus1');
            d1X2 = cwt(curve(sidx(end)).data(2,:),7,'gaus1');
            T = [d1X1;d1X2];
            N = [-d1X2;d1X1];
            sN = sum(N.*N,1).^-.5;
            sT = sum(T.*T,1).^-.5;
            T = bsxfun(@times,T,sT);
            N = bsxfun(@times,N,sN);
            initD = [N(:,tipIDX)';T(:,tipIDX)'];
            x = curve(sidx(end)).data(1,tipIDX);
            y = curve(sidx(end)).data(2,tipIDX);
            midline = trackFromPointAtGradient(sDIST,[x y]',initD(1,1),initD(1,2),1400,15,pi,[15 200],.3);
            
            
            
            %sample the euclidean distance along the midline
            sEDT = ba_interp2(sDIST,midline(1,:),midline(2,:));
            startSNIP = 15;
            endSNIP = 20;
            uE = mean(sEDT(startSNIP:startSNIP+endSNIP));
            sE = std(sEDT(startSNIP:startSNIP+endSNIP));
            kidx = find(sEDT > 15);
            midline = midline(:,1:kidx(1));
            
           
            if disp
                imshow(sI{i})
                hold on
                plot(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),'g','LineWidth',2);
                plot(midline(1,:),midline(2,:),'r');
                plot(curve(sidx(end)).data(1,tipIDX),curve(sidx(end)).data(2,tipIDX),'k*');
                drawnow
                pause(.1);
            end
            
            
            
            data(i).contourdata = sampleCurveBank(sI{i},curve(sidx(end)),NH,nhSZ);            
            data(i).midlinedata = midline;
            
            
            dist = [];
            for l = 1:size(data(i).contourdata.data,2)
                delta = bsxfun(@minus,data(i).contourdata.data(:,l),data(i).midlinedata);
                delta = sum(delta.*delta,1).^.5;
                dist(l) = min(delta);
            end
            fidx = find(dist < 15);            
            kidx = find(dist > 15);
            data(i).radical.data = data(i).contourdata.data(:,fidx);            
            data(i).kernel.data = data(i).contourdata.data(:,kidx);
            dL = diff(data(i).kernel.data,1,2);
            dL = sum(dL.*dL,1).^.5;
            [JUNK midx] = max(dL);
            data(i).kernel.data = circshift(data(i).kernel.data,[0 -midx]);
            midPoint = .5*(data(i).kernel.data(:,1) + data(i).kernel.data(:,end));
            data(i).kernel.data = [midPoint data(i).kernel.data midPoint];
            midPoint = .5*(data(i).radical.data(:,1) + data(i).radical.data(:,end));
            data(i).radical.data = [midPoint data(i).radical.data midPoint];
        end
        
        
        
        
        
        
        myQuickTracker(fileStack,CUR);
        
        p.type = 'disk';
        p.value{1} = [0 15 200];
        p.value{2} = [-pi pi 201];
        tD = phytoAdomain(p);
        tD.generateDomain();
        fun = @(x)generate_pbpca(x,fileStack,tD);                                
        af = distrib(CUR,fun);
        
        
        L = CUR.iLength();
        length = [];
        for i = 1:numel(L)
            length = [length L{i}];
        end
        
        if ~isempty(outPath)
            
            angle = af.angle();
            for e = 1:numel(angle)
                A(e,:) = angle{e};
            end
            [pth nm ext] = fileparts(fileName);
            fidx = strfind(pth,filesep);
            % logan spool
            csvwrite([outPath 'loganSPOOL' filesep pth(fidx(end)+1:end) '.csv'],A);
            
            % length spool
            LEVEL = 3;
            csvwrite([outPath 'csvSPOOL' filesep strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--ln.csv'],length');
            
            % growthrate spool
            LEVEL = 3;
            csvwrite([outPath 'csvSPOOL' filesep strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--gr.csv'],diff(length,1,1)');
            
            % csv spool
            LEVEL = 3;
            csvwrite([outPath 'csvSPOOL' filesep strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--an.csv'],A');
        end
        
        
        
        
        
        if disp
            close all
            h = imshow(IB,[]);
            hold on
            CUR.view(h,[]);
            af.view(h,[]);
            drawnow
            hold off
            [pth nm ext] = fileparts(fileName);
            fidx = strfind(pth,filesep);
            saveas(h,['/mnt/scratch5/maizeJUNK/' pth(fidx(end):end) '.tif']);
        end
        %}
        
    catch ME
        fprintf(['error@' fileName '\n']);
        error('test');
    end    
end