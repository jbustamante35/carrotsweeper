classdef maizeImage < handle
    properties
        filename;
        imageData;
    end
    
    methods
        function [obj] = maizeImage(filename)
            obj.filename = filename;
        end
        
        function [I] = imread(obj)
            obj.loadImage();
            I = obj.imageData;
        end
        
        function [] = loadImage(obj)
            if isempty(obj.imageData)
                obj.imageData = double(imread(obj.filename));
            end
        end

        function [kernelCenters] = extractKernelCenterPoints(obj)
            obj.loadImage();
            
            
            
            
            GLOBAL_UPPERLEFT = [10 300];
        

            I = double(obj.imageData(GLOBAL_UPPERLEFT(1):end-10,GLOBAL_UPPERLEFT(2):end-500))/255;
            % for non max
            fS2 = 31;
            fS1 = 201;
            fS = 31;                % smoothing for measure on finding (A)
            gradDiskSize = 11;
            CROPBOX_HALF_WIDTH = 100;

            % row samples width
            rowS = 20;              % thickness for sampling along rows and finding the cap
        
            
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
            kernelCenters = [cx';cy'];
            kernelCenters = bsxfun(@plus,kernelCenters,flipud(GLOBAL_UPPERLEFT'));
        end
        
        function [curveSet] = extractClosedContours(obj,containerClass,contourLevel,filterLevel,filterLength)
            if nargin == 4;filterLength=[0 inf];end
            obj.loadImage();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reduce image to gray scale if color
            if size(obj.imageData,3) == 3
                obj.imageData = rgb2gray(obj.imageData/255);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % smooth the image
            tmp = imfilter(obj.imageData,fspecial('gaussian',5*filterLevel,filterLevel),'replicate');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % obtain the gradient for contouring
            [g1 g2] = gradient(tmp);
            G = (g1.^2 + g2.^2).^.5;
            C = contourc(G,contourLevel);
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % select closed curves and filter on length
            [curve] = selectClosedCurves(curve);
            ridx = [curve.length] < filterLength(1) | [curve.length] > filterLength(2);
            curve(ridx) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % pour the curves into the given container class and assign the
            % file from which they came
            for e = 1:numel(curve)
                curveSet(e) = containerClass(curve(e).data);
                curveSet(e).filename = obj.filename;
            end
        end
        
        function [curveSet] = sampleGradient(obj,curves)
            obj.loadImage();
            [g1 g2] = gradient(obj.imageData);
            G = (g1.^2 + g2.^2).^.5;
            for e = 1:numel(curves)
                sample = ba_interp2(G,curves(e).segment(1,:),curves(e).segment(2,:));
                curves(e).totalGradient = mean(sample);
            end
        end
        
        function [I] = readCropped(obj,BOX)
            ROWS = round([BOX(1) BOX(1)+BOX(3)]);
            COLS = round([BOX(2) BOX(2)+BOX(4)]);
            imageData = double(imread(obj.filename,'PixelRegion',{COLS ROWS}));
            I = croppedMaizeImage( obj.filename);            
            I.cropBox = BOX;
            I.imageData = imageData;
        end
    end
end