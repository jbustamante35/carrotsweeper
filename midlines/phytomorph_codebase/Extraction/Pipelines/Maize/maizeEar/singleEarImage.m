function [KernelLength sM] = singleEarImage(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,addcut,baselineBlue,fill,CHUNK,windowSize,toSave,toDisplay)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                singleEarImage.m is main function to handle ear analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                StoN.m, checkBlue.m, measureKernelLength.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fileName:               An image to be analyze in a string that includes path and file name.
                noe:                    Number of cobs that are expected to be analyzed. 
                oPath:                  A path to result of analysis.
                rPath:                  iPlant location to save results.
                rawImage_scaleFactor:   A desired percentage to resize the image.
                checkBlue_scaleFactor:  A desired percentage to resiThe baseline threshold to remove blue header in checkBlue.
                fill:                   The radius of disk for Kernel of an image close operation.
                CHUNK:                  The number of chunk for input for FFT in myBlock0.
                windowSize:             The value for window size.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:              0 - not to save, 1 - toze the image in checkBlue.
                defaultAreaPix:         The default pixel to be considered noise relative to 1200 dpi.
                addcut:                 The boarder handle for checkBlue. This is an addition to blue top computed in checkBlue.
                baselineBlue:           The baseline threshold to remove blue header in checkBlue.
                fill:                   The radius of disk for Kernel of an image close operation.
                CHUNK:                  The number of chunk for input for FFT in myBlock0.
                windowSize:             The value for window size.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:              0 - not to save, 1 - to save.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    fprintf(['****************************************************************************************\n']);
    versionString = 'Publication Version 1.0 - Monday, March 28, 2016. \n';
    startString = 'Starting ear analysis algorithm. \n';
    fprintf([startString,versionString]);
    fprintf(['****************************************************************************************\n']);
    totalTimeInit = clock;
    %%%%%%%%%%%%%%%%%%%%%%%
    % init return vars    
    sM = [];   
    KernelLength = [];
    fprintf(['starting with variable and environment initialization.\n']);
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
        noe = StoN(noe);
        checkBlue_scaleFactor = StoN(checkBlue_scaleFactor);
        rawImage_scaleFactor = StoN(rawImage_scaleFactor);
        addcut = StoN(addcut);
        defaultAreaPix = StoN(defaultAreaPix);
        baselineBlue = StoN(baselineBlue);
        fill = StoN(fill);
        CHUNK = StoN(CHUNK);
        windowSize = StoN(windowSize);
        %%%%%%%%%%%%%%%%%%%%%%%
        % print out the fileName, number of ears, output path
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['FileName:' fileName '\n']);
        fprintf(['Number of Ears:' num2str(noe) '\n']);
        fprintf(['OutPath:' oPath '\n']); 
        fprintf(['Image resize in checkBlue:' num2str(checkBlue_scaleFactor) '\n']); 
        fprintf(['Raw image resize:' num2str(rawImage_scaleFactor) '\n']);  
        fprintf(['Threshold noise size:' num2str(defaultAreaPix) '\n']);  
        fprintf(['The window size:' num2str(windowSize) '\n']);
        fprintf(['The boarder handle for checkBlue:' num2str(addcut) '\n']);
        fprintf(['Baseline threshold to remove blue header:' num2str(baselineBlue) '\n']);
        fprintf(['The radius of disk for closing:' num2str(fill) '\n']);
        fprintf(['The number of chunk of blocks for FFT:' num2str(CHUNK) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
        [pth nm ext] = fileparts(fileName);
        fprintf(['ending with variable and environment initialization.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % read the image and take off the blue strip for bar codedefaultAreaPix,fracDpi
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%windowSizeD = StoN(windowSizeD);%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
        fprintf(['starting with image analysis. \n']);
        % the number of down sample grid sites
        gridSites = 10;
        [KernelLength FT BB S MT] = measureKernelLength(I,noe,windowSize,gridSites,defaultAreaPix,fill,CHUNK);
        % average kernel height
        uT = nanmean(KernelLength,2);        
        DATA = [];             
        % gather the data together for the bounding box
        for b = 1:numel(BB)            
            DATA = [DATA;[BB{b}(3:4) S.average_WIDTH(b)]];
        end
        % for new format for saving
        DATA = [DATA uT];
        DATA = reshape(DATA',[1 numel(DATA)]);
        fprintf(['ending with image analysis. \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
        if toDisplay
            fprintf(['starting with image and results display \n']);
            %%%%%%%%%%%%%%%%%%%%%%%
            % display the results
            %%%%%%%%%%%%%%%%%%%%%%%
            h = image(I);
            hold on
            %%%%%%%%%%%%%%%%%%%%%%%
            % make rectangle around the cob and plot the cosine function
            %%%%%%%%%%%%%%%%%%%%%%%
            for b = 1:numel(BB)
                rectangle('Position',BB{b},'EdgeColor','r');
                CS = [1:BB{b}(4)] + BB{b}(2);
                Func = 100*cos(2*pi/uT(b)*CS) + BB{b}(3)/2 + BB{b}(1);
                plot(Func,CS,'r');                
            end
            title(['Kernel Length:' num2str(nanmean(KernelLength,2)')]);
            %%%%%%%%%%%%%%%%%%%%%%%
            % format the return image
            %%%%%%%%%%%%%%%%%%%%%%%
            axis equal;axis off;drawnow;set(gca,'Position',[0 0 1 1]);
            fprintf(['ending with image and results display \n']);
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
            %% spin-up JSON output format
            %%%%%%%%%%%%%%%%%%%%%%%%%
            linkTo = stripiTicket(fileName);
            linkPath = stripiTicket(rPath);
            [jP,tN] = fileparts(fileName);
            for e = 1:numel(BB)
                
                tmpDoc = [];
                tmpDoc = generatePhenotypeNode(tmpDoc,linkTo,'orginalImage','orginalImage');
                tmpDoc = generatePhenotypeNode(tmpDoc,[linkPath filesep tN '_thumb.tif' ],'thumbNail','thumbNail');
                tmpDoc = generatePhenotypeNode(tmpDoc,S.widthProfile(e,:),{'along','position'},'widthProfile');
                %tmpDoc = generatePhenotypeNode(tmpDoc,S.RGB(e,:),{'along','colorIndex'},'rgb_color');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(3),{'maxWidth'},'maxWidth');
                tmpDoc = generatePhenotypeNode(tmpDoc,S.average_WIDTH(e),{'averageWidth'},'averageWidth');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(4),{'length'},'length');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e},{'boundingBox'},'boundingBox');
                tmpDoc = generatePhenotypeNode(tmpDoc,uT(e),{'kernelLength'},'kernelLength');
                phDoc(e) = tmpDoc;
                
                %phDoc(e).linkTo = linkTo;
                %phDoc(e).thumbImage = [linkPath filesep tN '_thumb.tif' ];
                %phDoc(e).widthProfile = S.widthProfile(e,:);
                %phDoc(e).rgb_color = S.RGB(e,:);
                %phDoc(e).maxWidth = BB{e}(3);
                %phDoc(e).averageWidth = S.average_WIDTH(e);
                %phDoc(e).length = BB{e}(4);
                %phDoc(e).boundingBox = BB{e};
                %phDoc(e).kernelLength = uT(e);
            end
            JSON_string = savejson('earDoc',phDoc);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save image file
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{1} = [oPath nm '_result.tif'];
            saveas(h,fileList{1});
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save mat file
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{2} = [oPath nm '.mat'];
            save(fileList{2},'BB','fileName','KernelLength','FT','S','windowSize');
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save csv data
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{3} = [oPath nm 'compile_results.csv'];
            csvwrite(fileList{3},DATA);
            fileList{4} = [oPath nm 'width_results.csv'];
            csvwrite(fileList{4},S.widthProfile);
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
        getReport(ME);
        fprintf(['******error in:singleEarImage.m******\n']);
    end
    close all
    fprintf(['****************************************************************************************\n']);
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
    endString = 'Ending ear analysis algorithm. ';
    fprintf([endString,versionString]);
    fprintf(['****************************************************************************************\n']);
end