function [S] = singleCobImage(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                singleCobImage.m is main function to handle cob analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box and color circle. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                StoN.m, checkBlue.m, maizeCob.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fileName:               An image to be analyze in a string that includes path and file name.
                noe:                    Number of cobs that are expected to be analyzed. 
                oPath:                  A path to result of analysis.
                rPath:                  iPlant location to save results.
                checkBlue_scaleFactor:  A desired percentage to resize the image in checkBlue.
                rawImage_scaleFactor:   A desired percentage to resize the image.
                defaultAreaPix:         The default pixel to be considered noise relative to 1200 dpi.
                rho:                    The radius of color circle, relative to 1200 dpi.
                addcut:                 The boarder handle for checkBlue. This is an addition to blue top computed in checkBlue.
                baselineBlue:           The baseline threshold to remove blue header in checkBlue.
                colRange1:              The color range for back ground to be removed in getcobMask.
                colRange2:              The color range for back ground to be removed in getcobMask.
                fill:                   The radius of disk for Kernel of an image close operation.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:      0 - not to save, 1 - to save.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

    algorithmName = 'maize-cobs';
    algorithmVerion = '1.0';
    fprintf(['****************************************************************************************\n']);
    versionString = 'Publication Version 1.0 - Monday, March 28, 2016. \n';
    startString = 'Starting cob analysis algorithm. \n';
    fprintf([startString,versionString]);
    fprintf(['****************************************************************************************\n']);
    totalTimeInit = clock;
    try 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with variable and environment initialization.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % init reporting for %uncLog
        %%%%%%%%%%%%%%%%%%%%%%%
        [rawName ticket] = stripiTicket(fileName);
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'2'},{'1'},1);
        %%%%%%%%%%%%%%%%%%%%%%%
        % init the icommands
        %%%%%%%%%%%%%%%%%%%%%%%
        %initIrods();
        %%uncLog({rawName},'set',algorithmName,algorithmVerion,{'3'},{'1'},1);
        %%%%%%%%%%%%%%%%%%%%%%%
        % convert the strings to numbers if they are strings
        %%%%%%%%%%%%%%%%%%%%%%%
        noe = StoN(noe);
        checkBlue_scaleFactor = StoN(checkBlue_scaleFactor);
        rawImage_scaleFactor = StoN(rawImage_scaleFactor);
        defaultAreaPix = StoN(defaultAreaPix);       
        rho = StoN(rho);
        addcut = StoN(addcut);    
        baselineBlue = StoN(baselineBlue);
        colRange1 = StoN(colRange1);
        colRange2 = StoN(colRange2);
        fill = StoN(fill);
        %%%%%%%%%%%%%%%%%%%%%%%
        % print out the input variables
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['FileName:' fileName '\n']);
        fprintf(['Number of Ears:' num2str(noe) '\n']);
        fprintf(['OutPath:' oPath '\n']);     
        fprintf(['Image resize in checkBlue:' num2str(checkBlue_scaleFactor) '\n']); 
        fprintf(['Raw image resize:' num2str(rawImage_scaleFactor) '\n']);  
        fprintf(['Threshold noise size:' num2str(defaultAreaPix) '\n']);
        fprintf(['The radius of color circle:' num2str(rho) '\n']);
        fprintf(['The boarder handle for checkBlue:' num2str(addcut) '\n']);
        fprintf(['Baseline threshold to remove blue header:' num2str(baselineBlue) '\n']);
        fprintf(['Background Color Range I:' num2str(colRange1) '\n']);
        fprintf(['Background Color Range II:' num2str(colRange2) '\n']);
        fprintf(['The radius of disk for closing:' num2str(fill) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
        [pth nm ext] = fileparts(fileName);
        fprintf(['ending with variable and environment initialization.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%
        % read the image and take off the blue strip for bar code
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with image load.\n']);
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'3'},{'1'},1);
        I = imread(fileName);
        % remove 4th pane for some images
        I = I(:,:,1:3);
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'4'},{'1'},1);
        thumb = imresize(I,.25);
        % rawImage_scaleFactor to lower 'DPI' effect, by fraction
        % If resize factor is 1, do not excecute imresize
        if rawImage_scaleFactor ~= 1;I = imresize(I,rawImage_scaleFactor);end
        % check blue header and remove
        I = checkBlue(I,checkBlue_scaleFactor,addcut,baselineBlue);
        fprintf(['ending with image load.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        % main measurement code call
        fprintf(['starting with image analysis \n']);
        % convert to single
        %I = single(I)/255;
        % run the analysis
        [BB S] = maizeCob(I,noe,defaultAreaPix,colRange1,colRange2,fill);
        % stack the results from the bounding box
        DATA = [];
        for b = 1:numel(BB)                
            DATA = [DATA;[BB{b}(3:4) S.average_WIDTH(b)]];        
        end
        uDATA = mean(DATA,1);
        sDATA = std(DATA,1,1);
        DATA = reshape(DATA',[1 numel(DATA)]);
        fprintf(['ending with image analysis \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toDisplay
            fprintf(['starting with image and results display \n']);
            h = image(I);
            hold on
            % make rectangle around the cob
            for b = 1:numel(BB)
                rectangle('Position',BB{b},'EdgeColor','r');
                UR = BB{b}(1:2);
                UR(1) = UR(1) + BB{b}(3);
                rectangle('Position',[UR rho rho],'EdgeColor','none','Curvature',[1 1],'FaceColor',S.RGB(b,:)/255);
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            % format the image
            axis equal;axis off;drawnow;set(gca,'Position',[0 0 1 1]);
            fprintf(['ending with image and results display \n']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toSave
            %% data to database
            for e = 1:numel(BB)
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:length:' num2str(e)]},{num2str(BB{b}(4))},1);
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:width:' num2str(e)]},{num2str(BB{b}(3))},1);
            end
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:average_length:' num2str(e)]},{num2str(uDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_length:' num2str(e)]},{num2str(sDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            
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
                tmpDoc = generatePhenotypeNode(tmpDoc,S.RGB(e,:),{'along','colorIndex'},'rgb_color');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(3),{'maxWidth'},'maxWidth');
                tmpDoc = generatePhenotypeNode(tmpDoc,S.average_WIDTH(e),{'averageWidth'},'averageWidth');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(4),{'length'},'length');
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e},{'boundingBox'},'boundingBox');
                phDoc(e) = tmpDoc;
            end
            JSON_string = savejson('cobDoc',phDoc);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% csv data file to disk
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList = {};
            fprintf(['starting save phase \n']);
            % save image
            fileList{end+1} = [oPath nm '_thumb.tif' ];
            imwrite(thumb,fileList{end});
            % save image
            fileList{end+1} = [oPath nm '_result.tif' ];
            saveas(h,fileList{end});
            % save mat file
            fileList{end+1} = [oPath nm '.mat'];
            save(fileList{end},'BB','fileName','S');
            % save the global parameters in file
            fileList{end+1} = [oPath nm '.csv'];
            csvwrite(fileList{end},DATA);
            % save the width parameters
            fileList{end+1} = [oPath nm '_width_results.csv'];
            csvwrite(fileList{end},S.widthProfile);
            % save the color values
            fileList{end+1} = [oPath nm '_cobRGB.csv'];
            csvwrite(fileList{end},S.RGB);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save JSON string
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm '_jdoc.json'];
            fileID = fopen(fileList{end},'w');
            fprintf(fileID,strrep(JSON_string,'\/','\\/'));
            fclose(fileID);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % push to iRODS
            %%%%%%%%%%%%%%%%%%%%%%%%%
            pushToiRods(rPath,fileList);
            fprintf(['ending save phase \n']);
        end
        close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'1'},1);
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:singleCobImage.m******\n']);
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'2'},1);
    end
    close all
    fprintf(['****************************************************************************************\n']);
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
    endString = 'Ending cob analysis algorithm. \n';
    fprintf([endString,versionString]);
    fprintf(['****************************************************************************************\n']);
    %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'6'},{'1'},1);
end

