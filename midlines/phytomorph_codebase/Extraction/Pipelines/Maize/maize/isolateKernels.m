function [] = isolateKernels(fileStack,disp,outPath)
    try
        warning off
        mkdir(outPath);
        
        GLOBAL_UPPERLEFT = [10 300];
        straight = 1;
        
        
        
        fileName = fileStack{1};
        I = imread(fileName);        
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
            
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find the highest curvature for the longest curve
            % first smooth large and find the biggest peak in the curvature
            % profile then find the curve within a window which is
            % contributing most the the bend
            LOC = bsxfun(@plus,curve(sidx(end)).data,CROPBOX{i}(1:2)');
            LOC = bsxfun(@minus,LOC,[cx(i);cy(i)]);
            KSNIP = 50;
            out = cwtK_closed(curve(sidx(end)).data',{60});
            out.K = out.K.*double(LOC(1,:) > 0);
            [J tipIDX] = min(out.K);            
            out = cwtK_closed(curve(sidx(end)).data',{10});
            out.K = out.K.*double(LOC(1,:) > 0);
            [J fine_tipIDX] = min(out.K(tipIDX-KSNIP:tipIDX+KSNIP));
            % find the curve point
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
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % the width of the belief or the momentum is a function of percent 
            % of cutoff - delta is the distance from straight
            % 0 will be used for now - but feedback could happen here -
            % notes
            wsigma = .3;
            maxWidth = 20*pi/180;
            width = @(delta)maxWidth*normpdf(delta,0,wsigma)/normpdf(0,0,wsigma);
            %widthTest = linspace(0,1,100);
            %plot(widthTest,width(widthTest));
            %%%%%%%%%%%%%%%%%%
            % function for radial belief
            k = 1.5;
            alpha = 1;
            scale = 10;
            radial = @(x)wblpdf(x*scale^-1,alpha,k);
            %radialTest = linspace(0,20,100);
            %plot(radialTest,radial(radialTest));
            %%%%%%%%%%%%%%%%%%
            % function for angular belief - the delta is the width and x is
            % theta
            angle = @(x,delta)normpdf(x,0,width(delta));
            %%%%%%%%%%%%%%%%%%
            % belief function
            func = @(rad,rho,delta)angle(rad,delta).*radial(rho);
            
            

            pointDensity = [15 300];
            RHO = 15;
            RAD = pi;
            %[rho rad] = ndgrid(linspace(0,RHO,pointDensity(1)),linspace(RAD,-RAD,pointDensity(2)));
            %CAR = [rho(:)'.*cos(rad(:)');rho(:)'.*sin(rad(:)')];
            %PHI = [rho(:)';rad(:)'];
            %TEST = func(PHI(2,:),PHI(1,:),0);
            %TEST = reshape(TEST,pointDensity);
            %figure;
            %mesh(reshape(CAR(1,:),pointDensity),reshape(CAR(2,:),pointDensity),TEST);




            d1X1 = cwt(curve(sidx(end)).data(1,:),7,'gaus1');
            d1X2 = cwt(curve(sidx(end)).data(2,:),7,'gaus1');
            T = [d1X1;d1X2];
            N = [-d1X2;d1X1];
            sN = sum(N.*N,1).^-.5;
            sT = sum(T.*T,1).^-.5;
            T = bsxfun(@times,T,sT);
            N = bsxfun(@times,N,sN);
            
            %imshow(BW,[]);
            %hold on
            %quiver(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),N(1,:),N(2,:))
            %quiver(curve(sidx(end)).data(1,tipIDX),curve(sidx(end)).data(2,tipIDX),N(1,tipIDX),N(2,tipIDX),'r')
            
            
            initD = [N(:,tipIDX)';T(:,tipIDX)'];
            x = curve(sidx(end)).data(1,tipIDX);
            y = curve(sidx(end)).data(2,tipIDX);

            T = goT();
            T.setWfunction(func);
            T.setNhoodRho(RHO);
            T.setNhoodRad(RAD);
            T.setNhoodDensity(pointDensity);
            T.generateH();
            T.setPosition([x;y]);
            T.setImage(sDIST);
            T.setDirection(initD);
            T.walk(300);
            T.reparameterizeCurve();            
            %alongRIDGE = ba_interp2(sDIST,T.position(1,:),T.position(2,:));
            %targetValue = mean(alongRIDGE(20:40));
            %targetValueError = std(alongRIDGE(20:40));
            %fidx = find((alongRIDGE-targetValue) > targetValueError);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            if disp
                imshow(sI{i})
                plot(curve(sidx(end)).data(1,:),curve(sidx(end)).data(2,:),'g','LineWidth',2);
                plot(T.position(1,:),T.position(2,:),'b');
                plot(curve(sidx(end)).data(1,tipIDX),curve(sidx(end)).data(2,tipIDX),'k*');
                drawnow
                pause(.1);
            end
            % shift for P = phytoApoint(allP);            tipcentered
            %cur = circshift(curve(sidx(end)).data(:,1:end-1),[0 -tipIDX]);
            %plot(cur(1,1),cur(2,1),'ko');
            %out = cwtK_closed(cur',{5});
            %curveBank{end+1} = cur;
            
        end
        
        
        
        
        
        
        myQuickTracker(fileStack,CUR);
        
        p.type = 'disk';
        p.value{1} = [0 15 200];
        p.value{2} = [-pi pi 201];
        tD = phytoAdomain(p);
        tD.generateDomain();
        fun = @(x)generate_pbpca(x,fileStack,tD);                                
        af = distrib(CUR,fun);

        
        %fun= @(x)iLe
        %af = distrib(CUR,iLength);
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
            tp1 = [outPath 'loganSPOOL' filesep 'angle' filesep];
            mkdir(tp1)
            csvwrite([tp1 pth(fidx(end)+1:end) '.csv'],A);
            % spool
            tp2 = [outPath 'loganSPOOL' filesep 'growthrate' filesep];
            mkdir(tp2)
            csvwrite([tp2 pth(fidx(end)+1:end) '.csv'],diff(length,1,1)');
            % spool 
            tp3 = [outPath 'loganSPOOL' filesep 'length' filesep];
            mkdir(tp3)
            csvwrite([tp3 pth(fidx(end)+1:end) '.csv'],length');
            
            
            
            
            tp3 = [outPath 'csvSPOOL' filesep];
            mkdir(tp3)
            % length spool
            LEVEL = 3;
            csvwrite([tp3 strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--ln.csv'],length');
            
            % growthrate spool
            LEVEL = 3;
            csvwrite([tp3 strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--gr.csv'],diff(length,1,1)');
            
            % csv spool
            LEVEL = 3;
            csvwrite([tp3 strrep(pth(fidx(end-LEVEL)+1:end),filesep,'----') '--an.csv'],A');
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
            %saveas(h,['/mnt/scratch5/maizeJUNK/' pth(fidx(end):end) '.tif']);
        end

        
    catch ME
        fprintf(['error@' fileName '\n']);
        rethrow(ME);
    end    
end