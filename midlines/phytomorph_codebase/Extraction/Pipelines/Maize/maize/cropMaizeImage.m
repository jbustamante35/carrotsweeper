function [] = cropMaizeImage(fileName,dataSetName,oPath,bulkSave,miniSave,matFile)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set internal display
        dispi = 0;
        extractCurve = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load mat file - if mat file given then hmm run
        if nargin == 6
            load(matFile);
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make output path
        mkdir(oPath)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert if is char
        if isa(fileName,'char')
            fileName = cellStringConvert(fileName);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert if is char
        if isa(bulkSave,'char')
            bulkSave = str2num(bulkSave);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert if is char
        if isa(miniSave,'char')
            miniSave = str2num(miniSave);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create the mat out file
        outFile = [oPath dataSetName];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % upper left corner for cropping out seedlings
        SEEDLING_WINDOW = [50 200 50 100];        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OP(s) for extracting curve information
        OP.levels = linspace(0,.75,100);
        OP.filter_area = [200 5000];
        OP.displayCurve = 0;
        OP.waveletPara = 11;
        OP.RAD = 30;
        OP.nRAD = 30;
        OP.THETA = [-pi pi];
        OP.nTHETA = 200;
        OP.segmentLength = 31;
        OP.savePath = oPath;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop window from file
        I = cropSeedlingsFromFile(fileName{1},SEEDLING_WINDOW);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate cropboxes
        CROPBOX = generateCropBox(I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save cropbox, numberofseedlines,number of frames
        numberSeedlings = numel(CROPBOX);
        numberFrames = numel(fileName);
        if bulkSave | miniSave 
            save(outFile,'numberSeedlings','numberFrames','CROPBOX');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make movies
        for s = 1:numberSeedlings
            % init vars
            close all
            tipAngle = [];
            curveLabels = [];
            tipIndex = [];
            tipVec = [];
            try
                for tm = 1:numberFrames
                        tim = clock;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % make crop window
                        WIN = {[CROPBOX{s}(2)+SEEDLING_WINDOW(1) CROPBOX{s}(2)+CROPBOX{s}(4)+SEEDLING_WINDOW(1)]...
                               [CROPBOX{s}(1)+SEEDLING_WINDOW(2) CROPBOX{s}(1)+CROPBOX{s}(3)+SEEDLING_WINDOW(2)]};
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % read in frame
                        frame = double(imread(fileName{tm},'PixelRegion',WIN))/255;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % bindvector
                        frame = bindVec(frame);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if extractCurve
                            % get contour
                            curve = collectContourInformation_ver3(frame,OP);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if loaded hmm and pcaData then labelRoots
                        if nargin == 6
                            skelPoints = [];
                            % label contour and get tip angle
                            [curveLabels{tm} tipIndex(tm) tipVec(:,tm) tipAngle(tm)] = labelRoots(curve,pcaData,hmm,@(x)min(x),4,2,30);
                            % added to trace out midline
                            fidx = find(curveLabels{tm}==2 | curveLabels{tm}==6);
                            basePoint = mean(curve(1).data(:,fidx),2);
                            MASK = poly2mask(curve(1).data(1,:),curve(1).data(2,:),size(frame,1),size(frame,2));
                            skel = bwmorph(MASK,'skel',inf);
                            [skelPoints(:,1),skelPoints(:,2)] = find(skel);
                            startPoint = curve(1).data(:,tipIndex(tm));
                            path = traceMaizeMidline(skelPoints,startPoint',basePoint');
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if bulkSave
                            % save data
                            saveBulkData(tm,s,frame,curve,outFile);
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if miniSave
                            curve.S = [];
                            curve.E = [];
                            curve.segs = [];
                            % create curve string
                            curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                            cmd = [curveString '=curve;'];                
                            eval(cmd);
                            save(outFile,curveString,'-append');
                            % create path string
                            pathString = ['seedling' num2str(s) 'midline' num2str(tm)];
                            cmd = [pathString '=path;'];
                            eval(cmd);
                            save(outFile,pathString,'-append');
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % if display
                        if dispi
                            CL = {'r' 'b' 'g' 'm' 'c' 'y' 'r' 'b' 'g' 'm' 'c' 'y'};
                            % show results
                            imshow(frame,[]);
                            hold on
                            plot(curve(1).data(1,:),curve(1).data(2,:),'r');
                            UQg = unique(curveLabels{tm});
                            for g = 1:numel(UQg)
                                fidx = find(curveLabels{tm}==UQg(g));
                                plot(curve(1).data(1,fidx),curve(1).data(2,fidx),[CL{UQg(g)} '.'],'MarkerSize',1);
                            end
                            plot(curve(1).data(1,tipIndex(tm)),curve(1).data(2,tipIndex(tm)),'k*');
                            plot(basePoint(1),basePoint(2),'y*')
                            plot(path(2,:),path(1,:),'g')
                            hold off
                            drawnow
                            
                        end
                        fprintf(['Complete with seedling:frame@' num2str(s) ':' num2str(tm) '@' num2str(etime(clock,tim)) '\n']);

                end
                if miniSave
                    % save data
                    saveMiniData(s,curveLabels,tipAngle,tipIndex,tipVec,outFile);
                end
            catch ME
                % report error
                disp(ME.getReport());
                % init vars to zero
                curveLabels{1} = [];
                tipAngle = zeros(1,numberFrames);
                tipIndex = zeros(1,numberFrames);
                tipVec = zeros(2,numberFrames);
                if miniSave
                    % save data
                    saveMiniData(s,curveLabels,tipAngle,tipIndex,tipVec,outFile);
                end
            end
        end
    catch ME
        disp(ME.getReport);
    end
      
end



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
    % added to remove the case when there is only one kernel
    % and the kernel does have a stronger signal than the edge
    fidx(fidx > .75*size(I,2)) = [];
    % added above
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

function [I] = cropSeedlingsFromFile(fileName,SEEDLING_WINDOW)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read in the window which should contain the seedlings
        I = imread(fileName);
        I = double(I(SEEDLING_WINDOW(1):end-SEEDLING_WINDOW(3),SEEDLING_WINDOW(2):end-SEEDLING_WINDOW(4)))/255;
end

function [] = saveBulkData(tm,s,frame,curve,outFile)
    fprintf(['Spoolding seedling:frame@' num2str(s) ':' num2str(tm) '\n'])
    frameString = ['seedling' num2str(s) 'frame' num2str(tm)];                
    cmd =  [frameString '=frame;'];
    eval(cmd);
    curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
    cmd = [curveString '=curve;'];                
    eval(cmd);
    if ~exist(outFile)
        save(outFile,frameString,curveString);
    else
        save(outFile,frameString,curveString,'-append');
    end
end

function [] = saveMiniData(s,curveLabels,tipAngle,tipIndex,tipVec,outFile)
    fprintf(['Spoolding mini seedling@' num2str(s) '\n'])
    % create curve label string
    curveLabelString = ['seedling' num2str(s) 'curveLabels'];
    cmd = [curveLabelString '=curveLabels;'];                
    eval(cmd);
    % create angle string
    angleString = ['seedling' num2str(s) 'angle'];
    cmd = [angleString '=tipAngle;'];                
    eval(cmd);
    % create index string
    tipIndexString = ['seedling' num2str(s) 'tipIndex'];
    cmd = [tipIndexString '=tipIndex;'];                
    eval(cmd);
    % create angle string
    tipVecString = ['seedling' num2str(s) 'tipVec'];
    cmd = [tipVecString '=tipVec;'];                
    eval(cmd);
    
    if ~exist(outFile)
        save(outFile,curveLabelString,angleString,tipIndexString,tipIndexString,tipVecString);
    else
        save(outFile,curveLabelString,angleString,tipIndexString,tipIndexString,tipVecString,'-append');
    end
end