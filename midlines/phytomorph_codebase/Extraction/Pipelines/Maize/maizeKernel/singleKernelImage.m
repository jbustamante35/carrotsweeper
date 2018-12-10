function [] = singleKernelImage(fileName,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,addcut,baselineBlue,fill,toSave,toDisplay)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                singleKernelImage.m is main function to handle kernel analysis. It takes all input variables 
                for its dependent functions. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                StoN.m, checkBlue.m, getThresholdLevel.m, getInitialGuessForTip.m,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fileName:               An image to be analyze in a string that includes path and file name.
                oPath:                  A path to result of analysis.
                rPath:                  iPlant location to save results.
                rawImage_scaleFactor:   A desired percentage to resize the image.
                checkBlue_scaleFactor:  A desired percentage to resize the image in checkBlue.
                defaultAreaPix:         The default pixel to be considered noise relative to 1200 dpi.
                addcut:                 The boarder handle for checkBlue. This is an addition to blue top computed in checkBlue.
                baselineBlue:           The baseline threshold to remove blue header in checkBlue.
                fill:                   The radius of disk for Kernel of an image close operation.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:              0 - not to save, 1 - to save.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    fprintf(['****************************************************************************************\n']);
    versionString = 'Publication Version 1.0 - Monday, March 28, 2016. \n';
    startString = 'Starting kernel analysis algorithm. \n';
    fprintf([startString,versionString]);
    fprintf(['****************************************************************************************\n']);
    totalTimeInit = clock;
    %%%%%%%%%%%%%%%%%%%%%%%
    % init vars
    MAJOR = [];
    MINOR = [];
    toSaveContour = [];
    %%%%%%%%%%%%%%%%%%%%%%%
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with variable and environment initialization.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % init the icommands
        %%%%%%%%%%%%%%%%%%%%%%%
        %initIrods();
        %%%%%%%%%%%%%%%%%%%%%%%
        % convert the strings to numbers if they are strings
        %%%%%%%%%%%%%%%%%%%%%%%
        checkBlue_scaleFactor = StoN(checkBlue_scaleFactor);
        rawImage_scaleFactor = StoN(rawImage_scaleFactor);
        defaultAreaPix = StoN(defaultAreaPix);
        addcut = StoN(addcut);
        baselineBlue = StoN(baselineBlue);
        fill = StoN(fill);
        %%%%%%%%%%%%%%%%%%%%%%%
        % print out the fileName, number of ears, output path
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['FileName:' fileName '\n']);
        fprintf(['OutPath:' oPath '\n']);
        fprintf(['Raw image resize:' num2str(rawImage_scaleFactor) '\n']);
        fprintf(['Image resize in checkBlue:' num2str(checkBlue_scaleFactor) '\n']); 
        fprintf(['Threshold noise size:' num2str(defaultAreaPix) '\n']);
        fprintf(['The boarder handle for checkBlue:' num2str(addcut) '\n']);
        fprintf(['Baseline threshold to remove blue header:' num2str(baselineBlue) '\n']);
        fprintf(['The radius of disk for closing:' num2str(fill) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % get the file parts
        %%%%%%%%%%%%%%%%%%%%%%%
        [pth nm ext] = fileparts(fileName);
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
        fprintf(['starting with variable and environment initialization.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % read the image and take off the blue strip for bar code
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with image load.\n']); 
        I = imread(fileName);
        I = I(:,:,1:3);
        % rawImage_scaleFactor to lower 'DPI' effect, by fraction 
        % If resize factor is 1, do not excecute imresize
        if rawImage_scaleFactor ~= 1;I = imresize(I,rawImage_scaleFactor);end
        I = checkBlue(I,checkBlue_scaleFactor,addcut,baselineBlue);
        fprintf(['ending with image load.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
        fprintf(['starting with image analysis\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % main part of code - threshold,filter and count
        % get the threshold level
        [level,G] = getThresholdLevel(I);        
        % threshold the image
        B = single(G) > level;        
        % remove small objects
        B = bwareaopen(B,defaultAreaPix);                
        % fill holes
        B = imfill(B,8,'holes');
        % remove objects connected to the kernel and are thin
        B = imopen(B,strel('disk',fill,8));
        % fill holes
        B = imfill(B,4,'holes');
        % re-remove small objects
        B = bwareaopen(B,defaultAreaPix);        
        % measure region props
        R = regionprops(B,'Area','MajorAxis','MinorAxis','Image','Centroid','Orientation','Eccentricity','PixelIdxList');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove QR code
        % added Feb 27,2018
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Area = [R.Area];
        ridx = find(Area > 1500000);
        if ~isempty(ridx)
            B(R(ridx).PixelIdxList) = 0;
        end
        % remeasure region props
        R = regionprops(B,'Area','MajorAxis','MinorAxis','Image','Centroid','Orientation','Eccentricity','PixelIdxList');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % stack the area
        AREA = [R.Area];
        %[fidx sel uA] = count(AREA);
        [fidx] = count(AREA);
        AREA = AREA(fidx==1);
        % get the kernel count
        KC = sum(fidx);
        %%%%%%%%%%%%%%%%%%%%%%%
        % boundary analysis        
        dB = bwboundaries(B);
        % select boundaries of single kernels
        dB = dB(fidx==1);
        [tipPoint dB] = getInitialGuessForTip(dB);
        [dB] = findTipPoints_fast(dB,B,I);
        [M] = getKernelMeasurements(B,dB);
        MAJOR = [M.MajorLength];
        MINOR = [M.MinorLength];
        toSaveContour = [];
        for e = 1:numel(M)
            toSaveContour = [toSaveContour M(e).iContour(:)];
        end
        fprintf(['ending with image analysis\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:numel(M)
            BOX = point2Box(fliplr(M(e).kernelMean),[M(e).MinorLength*2 M(e).MinorLength*2]);
            H = imcrop(I,BOX);
            dSeg = diff(M(e).MinorLine,1,1);
            A = atan2(-dSeg(2),dSeg(1));
            RM = [[cos(A) -sin(A)];[sin(A) cos(A)]];
            rC = (M(e).iContour);
            rI = imrotate(H,-A*180/pi);
           
            
            sz = size(rI);
            sz(3) = [];
            sz = fliplr(sz)/2;
            rC = bsxfun(@plus,rC,sz);
            
            
            
            MSK = poly2mask(rC(:,2),rC(:,1),size(rI,2),size(rI,1));
            R = regionprops(logical(MSK),'BoundingBox','Centroid');
            DELTA = R(1).Centroid - [size(rI,1) size(rI,2)]/2;
            GG = mean(bbox2points(R(1).BoundingBox));
            DELTA = R(1).Centroid - GG;
            R(1).BoundingBox(1:2) = R(1).BoundingBox(1:2) - 10;
            R(1).BoundingBox(3:4) = R(1).BoundingBox(3:4) + 20;
            rI = imcrop(rI,R(1).BoundingBox);
            
            
            %{
            imshow(rI,[]);
            rC = (M(e).iContour);
            sz = size(rI)/2;
            sz(3) = [];
            %sz = fliplr(sz);
            rC = bsxfun(@plus,rC,sz);
            rC = bsxfun(@plus,rC,flipdim(DELTA,2));
            %}
            ROT{e} = rI;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toDisplay
            renderResults(M,I,KC);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if toSave
            fprintf(['starting save phase \n']);
            
            
            
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save image file
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{1} = [oPath nm '.tif'];
            saveas(gca,fileList{1});
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save csv data
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{2} = [oPath nm '-individuals.csv'];
            csvwrite(fileList{2},[MAJOR' MINOR' AREA']);
            fileList{3} = [oPath nm '.csv'];
            csvwrite(fileList{3},[mean(MAJOR) std(MAJOR) mean(MINOR) std(MINOR) mean(AREA) std(AREA) KC]);
            fileList{4} = [oPath nm '-contours.csv'];
            csvwrite(fileList{4},toSaveContour');
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% spin-up JSON output format
            %%%%%%%%%%%%%%%%%%%%%%%%%
            linkTo = stripiTicket(fileName);
            linkPath = stripiTicket(rPath);
            [jP,tN] = fileparts(fileName);
            for e = 1:size(toSaveContour,2)
                
                tmpDoc = [];
                tmpDoc = generatePhenotypeNode(tmpDoc,linkTo,'orginalImage','orginalImage');
                tmpDoc = generatePhenotypeNode(tmpDoc,[linkPath filesep tN '_thumb.tif' ],'thumbNail','thumbNail');
                tmpDoc = generatePhenotypeNode(tmpDoc,M(e).kernelReferenceFrame,{'dim','basisVector'},'kernelReferenceFrame');
                tmpDoc = generatePhenotypeNode(tmpDoc,M(e).kernelMean,{'centroid','dim'},'kernelCentroid');
                tmpDoc = generatePhenotypeNode(tmpDoc,M(e).iContour,{'along','dim'},'contour');
                tmpDoc = generatePhenotypeNode(tmpDoc, MAJOR(e),{'majorAxis'},'majorAxis');
                tmpDoc = generatePhenotypeNode(tmpDoc,MINOR(e),{'minorAxis'},'minorAxis');
                tmpDoc = generatePhenotypeNode(tmpDoc,AREA(e),{'area'},'kernelArea');
                phDoc(e) = tmpDoc;
                
                %{
                phDoc(e).contourRasterFormat = toSaveContour(:,e)';
                phDoc(e).referenceFrame = M(e).kernelReferenceFrame;
                phDoc(e).kernelCentroid = M(e).kernelMean;
                phDoc(e).XYcontour = M(e).iContour;
                phDoc(e).majorAxis = MAJOR(e);
                phDoc(e).minorAxis = MINOR(e);
                phDoc(e).area = AREA(e);
                %}
            end
            JSON_string = savejson('kernelDoc',phDoc);
            %%%%%%%%%%%%%%%%%%%%%%%%%
               
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save JSON string
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm '_jdoc.json'];
            fileID = fopen(fileList{end},'w');
            fprintf(fileID,strrep(JSON_string,'\/','\\/'));
            fclose(fileID);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save image file
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm 'kernelModelData.mat'];
            save(fileList{end},'ROT');
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % push to iRODS
            %%%%%%%%%%%%%%%%%%%%%%%%%
            pushToiRods(rPath,fileList);
            fprintf(['ending save phase \n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%
        end
        close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:singleKernelImage.m******\n']);
    end
    close all
    fprintf(['****************************************************************************************\n']);
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
    endString = 'Endinging kernel analysis algorithm. \n';
    fprintf([endString,versionString]);
    fprintf(['****************************************************************************************\n']);
end