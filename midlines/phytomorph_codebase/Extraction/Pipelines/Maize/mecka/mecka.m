function [] = mecka(algorithm,fileName,numberOfObjects,oPath,remotePath,toSave,toDisplay,scanResolution,rawImage_scaleFactor)
    % warnings will not be printed on screen.
    warning off
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                mecka.m is wrapper function to handle Maize Ear Cob Kernel Analysis. A user is 
                to choose one of three algorithm to use. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                StoN.m, singleCobImage.m, singleEarImage.m, singleKernelImage.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                algorithm:              The argorithm to use. 'c' for singleCobImage.m, 'e' for singleEarImage.m, and 'k' for singleKernelImage.m
                fileName:               An image to be analyze in a string that includes path and file name.
                numberOfObjects:        Number of cobs that are expected to be analyzed. 
                oPath:                  A path to result of analysis.
                remotePath:             iPlant location to save results.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:              0 - not to save, 1 - to save.
                scanResolution:         The resolution that image was scanned in DPI. 
                rawImage_scaleFactor:   A desired percentage to resize the image. If no input, default is 1, noresizing.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        if ~strcmp(oPath(end),filesep)
            oPath = [oPath filesep];
        end
        % no input of rawImage_scaleFactor is default
        if nargin ~= 8
            % set to default value of 1
            rawImage_scaleFactor = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% common arguments for Ear, Cob and Kernel %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set to default value of 1200
        defaultResolution = 1200;
        scanResolution = StoN(scanResolution);
        % compute proportion of resolution over default
        fracDpi = scanResolution/defaultResolution;
        % set to default value of .25
        checkBlue_scaleFactor = .25; 
        % set to default value of 100
        addcut = 150;
        addcut = round(addcut*fracDpi);
        % set to default value of 600
        baselineBlue = 600;
        baselineBlue = round(baselineBlue*fracDpi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch algorithm
            case 'e'
                % set to default value of 10
                CHUNK = 10;
                % set to default value of 10^6
                defaultAreaPix = 10^6;
                defaultAreaPix = round(defaultAreaPix*fracDpi);
                % set to default value of 31
                fill = 31;
                % set to default value of 1200:25:1600
                windowSize = round(1200*fracDpi):round(25*fracDpi):round(1600*fracDpi);
                [KernelLength sM] = singleEarImage(fileName,numberOfObjects,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,addcut,baselineBlue,fill,CHUNK,windowSize,toSave,toDisplay);
            case 'c'
                % set to default value of 10^6
                defaultAreaPix = 10^6;
                defaultAreaPix = round(defaultAreaPix*fracDpi);
                % set to default value of 300
                rho = 300;
                rho = round(rho*fracDpi);
                % set to default value of 70
                colRange1 = 70;
                % set to default value of 166
                colRange2 = 166;
                % set to default value of 50
                fill = 50;
                singleCobImage(fileName,numberOfObjects,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay);
            case 'k'
                % numberOfObjects is not in use for kernel
                numberOfObjects;
                % set to default value of 50000
                %defaultAreaPix = 50000;
                defaultAreaPix = 5000; % changed for mockerlab 
                defaultAreaPix = round(defaultAreaPix*fracDpi);
                % set to default value of 31
                fill = 31;
                singleKernelImage(fileName,oPath,remotePath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,addcut,baselineBlue,fill,toSave,toDisplay);                
        end
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:mecka.m******\n']);
    end
end