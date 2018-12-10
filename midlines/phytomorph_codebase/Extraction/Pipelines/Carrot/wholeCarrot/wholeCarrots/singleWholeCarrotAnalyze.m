function [] = singleWholeCarrotAnalyze(filename,chop,crownLength,PETthresh,oPath,rPath)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                singleWholeCarrotAnalyze.m is main function to handle carrot image analysis.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                unsure as of 02.18.2017
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fileName:               An image to be analyze in a string that includes path and file name.
                chop:                   The pixel value to chop the root from the shoot
                crownLength:
                oPath:                  The local path to write data to before pushing to CyVerse
                rPath:                  The remote pach on CyVerse to resulting data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    fprintf(['*************************************************************************************************\n']);
    versionString = ['Starting carrot analysis algorithm. \nPublication Version 1.0 - Monday, April 20, 2018. \n'];
    fprintf(versionString);
    fprintf(['*************************************************************************************************\n']);
    totalTimeInit = clock;
    
    
    try
        csvHEADER = {};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['*********************************************************************\n']);
        fprintf(['starting with variable and environment initialization.\n']);
        fprintf(['*********************************************************************\n']);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % convert the strings to numbers if they are strings
        %%%%%%%%%%%%%%%%%%%%%%%
        [chop] = StoN(chop);
        [crownLength] = StoN(crownLength);
        filePushList = {};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the result files names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
        [p,n,e] = fileparts(filename);
        image_top_oName = [oPath n '_top.tif'];
        % store for push to irods
        filePushList{end+1} = image_top_oName;
        image_root_oName = [oPath n '_root.tif'];
        % store for push to irods
        filePushList{end+1} = image_root_oName;
        % store the tech top file name for push to irods
        image_topTech_oName = [oPath n '_topTech.tif'];
        % store for push to irods
        filePushList{end+1} = image_topTech_oName;
        fprintf(['*********************************************************************\n']);
        fprintf(['ending with variable and environment initialization.\n']);
        fprintf(['*********************************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the result files names
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image, make mask, and flip to proper direction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['*********************************************************************\n']);
        fprintf(['starting mask generation.\n']);
        fprintf(['*********************************************************************\n']);
        EXPECTED_FILE_SIZE = [6016 4016 3];
        % read the image
        fprintf(['reading image.\n']);
        I = imread(filename);
        % make the root mask
        fprintf(['making mask.\n']);
        [MASK] = makeRootMask(I);
        % flip and rotate the image
        [I,MASK] = flipAndRotate(I,MASK,EXPECTED_FILE_SIZE);
        % select the largest object from the binary mask
        MASK = bwlarge(MASK);
        % chop the carrot
        [root,top] = chopCarrot(MASK,chop,crownLength);
        % remove small holes ~ imfill small holes
        topB = ~bwareaopen(~top,300);
        fprintf(['*********************************************************************\n']);
        fprintf(['ending mask generation.\n']);
        fprintf(['*********************************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image, make mask, and flip to proper direction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
       

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MEASURE BLOCK
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['*********************************************************************\n']);
        fprintf(['starting measure block.\n']);
        fprintf(['*********************************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get classic measurements
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['***********************************\n']);
        fprintf(['starting classic measurements.\n']);
        fprintf(['***********************************\n']);
        topB = bwlarge(topB);
        classic{1} = regionprops(topB,'ConvexHull','EquivDiameter','Eccentricity','EulerNumber','BoundingBox','Perimeter','Solidity');
        classic{2} = regionprops(root,'ConvexHull','EquivDiameter','Eccentricity','EulerNumber','BoundingBox','Perimeter','Solidity');
        classic{3} = regionprops(MASK,'ConvexHull','EquivDiameter','Eccentricity','EulerNumber','BoundingBox','Perimeter','Solidity');
        fprintf(['***********************************\n']);
        fprintf(['ending classic measurements.\n']);
        fprintf(['***********************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get classic measurements
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the top and bottom
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['***********************************\n']);
        fprintf(['starting total digital biomass measurements.\n']);
        fprintf(['***********************************\n']);
        % measure digial bio mass
        dbm = [sum(top(:)) sum(root(:))];
        fprintf(['***********************************\n']);
        fprintf(['ending total digital biomass measurements.\n']);
        fprintf(['***********************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the top and bottom
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the bottom
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['***********************************\n']);
        fprintf(['starting root specific.\n']);
        fprintf(['***********************************\n']);
        % root profile is 3200 numbers long and is zero padded to fill
        TOTROOT = 3200;
        [profile,profileTrace] = getRootProfile(root,image_root_oName);
        profile = [profile;zeros(3200-numel(profile),1)];
        fprintf(['***********************************\n']);
        fprintf(['ending root specific.\n']);
        fprintf(['***********************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the bottom
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the top
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['***********************************\n']);
        fprintf(['starting shoot specific measurements.\n']);
        fprintf(['***********************************\n']);
        % perform the top sweep for top measurements
        [sweepT,integrationSweep,petNUM,petDIA,offset] = topSweep(double(top),1,image_top_oName);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % petiole widths
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get distance transform
        dist = bwdist(~topB);
        % get skeleton
        skel = bwmorph(topB,'skeleton',inf);
        % remove spurs
        skel = bwmorph(skel,'spur',100);
        % re-skelton
        skel = bwmorph(skel,'skeleton',inf);
        % find the skeleton locations
        sidx = find(skel);
        % sample the distance transform on the midline
        sampDist = dist(sidx);
        % less than PETthresh should be petioles - filter
        kidx = find(sampDist < PETthresh);
        % make blank mask for petioles
        Z = zeros(size(skel));
        % fill in the selected areas only
        Z(sidx(kidx)) = 1;
        % remove small none petiole fork residules
        Z = bwareaopen(Z,30);
        % distance transform by skeleton
        G = double(Z.*dist);
        % integrate for pet widths
        [~,integrationSweep_PETONLY,~,~] = topSweep(G,0,[],offset);
        % find the pixels for the petiles
        petIdx = find(Z);
        % sample the distance transform at the pet locations
        petWidth = dist(petIdx);
        % get the average pet width
        petWidth = mean(petWidth);
        % integrate the pet only sweep
        petCountFunc = integrationSweep_PETONLY.*petWidth.^-1;
        % take the average along the radius
        newPetCount = mean(petCountFunc);
        % make count stat along the radius
        petCountFuncMask = petCountFunc*newPetCount^-1;
        SH = histcounts(round(sampDist),0:300);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % petiole widths
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % count number of which are greater than fraction of longest
        % measure the objects in the pet mask
        R = regionprops(Z,'Area','PixelIdxList');
        % the area approximate is the length
        LENish = [R.Area];
        % normalize to the longest(area) path
        LENish = LENish / max(LENish);
        % find those objects which are greater than 20 percent of longest path
        fidx = find(LENish > .2);
        % another measure of pet count is the number of objects selected
        petCount = numel(fidx);
        
        fprintf(['***********************************\n']);
        fprintf(['ending shoot specific measurements\n']);
        fprintf(['***********************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the top
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        fprintf(['*********************************************************************\n']);
        fprintf(['ending measure block.\n']);
        fprintf(['*********************************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MEASURE BLOCK
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISPLAY and MEASURE CODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        image(255*top);
        colormap(gray);
        hold on
        % display selected paths in yellow
        for e = 1:numel(fidx)
            [p1 p2] = ind2sub(size(skel),R(fidx(e)).PixelIdxList);
            plot(p2,p1,'y.');
        end
        
        BOTCEN = [size(skel,1) mean(find(top(end,:)))];
        % sample top THREE largest
        [J,sidx] = sort([R.Area],'descend');
        TOPN = 3;
        TOPN = min(numel(R),TOPN);
        for e = 1:TOPN
            petW(e) = mean(dist(R(sidx(e)).PixelIdxList));
            petL(e) = J(e);
            [p1 p2] = ind2sub(size(skel),R(sidx(e)).PixelIdxList);
            plot(p2,p1,'g.');
            [J1,midx] = max(p1);
            plot([p2(midx) BOTCEN(2)],[p1(midx) BOTCEN(1)],'b');
            plot(BOTCEN(2),BOTCEN(1),'b*');
            extra = (diff([p2(midx) BOTCEN(2)]).^2 + diff([p1(midx) BOTCEN(1)]).^2).^.5;
            petL(e) = petL(e) + extra;
        end
        petW = mean(petW);
        
        saveas(gca,image_topTech_oName);
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISPLAY and MEASURE CODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the header information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        csvHEADER = {};
        for e = 1:TOTROOT
            csvHEADER{end+1} = ['rootProfile_' num2str(e)];
        end
        for e = 1:numel(integrationSweep)
            csvHEADER{end+1} = ['integrationSweep_' num2str(e)];
        end
        csvHEADER{end+1} = 'shootDigitalBioMass';
        csvHEADER{end+1} = 'rootDigitalBioMass';
        for e = 1:numel(SH)
            csvHEADER{end+1} = ['histogram_' num2str(e)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the header information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % csv file format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % write out "classic" binary morphology data
        ClassicNames = {'top','bottom','whole'};
        f = fields(classic{1});
        for e = 1:3
            for fl = 1:numel(f)
                 csvOutName = [oPath n '_' ClassicNames{e} '_' f{fl} '.csv'];
                 filePushList{end+1} = csvOutName;
                 if ~isempty(classic{e}.(f{fl}))
                    csvwrite(csvOutName,classic{e}.(f{fl}));
                 end
            end
        end
        
        % write out root integration sweep
        csvOutName = [oPath n '_rootProfile.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,profile');
        
        % write out integration sweep file
        csvOutName = [oPath n '_topProfile.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,integrationSweep);
        
        % write shoot digital biomass
        csvOutName = [oPath n '_shootDigitalBioMass.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,dbm(1));
        
        % write root digital biomass
        csvOutName = [oPath n '_rootDigitalBioMass.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,dbm(2));
        
        % write out raw distance transform values on skeleton
        csvOutName = [oPath n '_rawSkeletonDistanceTransform.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,sampDist);
        
        % write out histogram for distance transform values on skeleton
        csvOutName = [oPath n '_histogramDistanceTransform.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,SH);
        
        % make mastervec
        masterVec = [profile' integrationSweep dbm SH];
        for e = 1:numel(masterVec)
            MVO{1,e} = csvHEADER{e};
            MVO{2,e} = masterVec(e);
        end
        csvOutName = [oPath n '_mastervec.csv'];
        filePushList{end+1} = csvOutName;
        cell2csv(csvOutName,MVO);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % csv file format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % quality control image(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out top mask
        csvOutName = [oPath n '_shootMask.tif'];
        filePushList{end+1} = csvOutName;
        imwrite(topB,csvOutName); 
        
        % write out bottom mask
        csvOutName = [oPath n '_rootMask.tif'];
        filePushList{end+1} = csvOutName;
        imwrite(root,csvOutName); 
        
        % write out whole mask
        csvOutName = [oPath n '_wholeMask.tif'];
        filePushList{end+1} = csvOutName;
        imwrite(MASK,csvOutName);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % quality control image(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % strip ticket
        linkTo = stripiTicket(filename);
        linkPath = stripiTicket(rPath);
        [jP,tN] = fileparts(filename);
        % generate json style format document
        phenoTypeDocument = [];
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,linkTo,'orginalImage','orginalImage');
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,profile,{'along','position'},'rootWidthProfile');
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,integrationSweep,{'along','position'},'shootDigialBiomassProfile');
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,dbm(1),{'topTotalBiomass'},'topTotalDigitalBiomass');
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,dbm(2),{'bottomTotalBiomass'},'rootTotalDigitalBiomass');
        phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,SH,{'along','position'},'distanceTransformHistogram');
        for e = 1:numel(ClassicNames)
            flds = fields(classic{e});
            flds = setdiff(flds,{'ConvexHull'});
            for f = 1:numel(flds)
                phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,classic{e}.(flds{f}),flds{f},flds{f});
            end
        end
        % generate json string
        JSON_string = savejson('carrotDoc',phenoTypeDocument);
        % save json document
        filePushList{end+1} = [oPath n '_jdoc.json'];
        fileID = fopen(filePushList{end},'w');
        fprintf(fileID,strrep(JSON_string,'\/','\\/'));
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % old JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate json string - OLD format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phDoc.rootProfile = profile';
        phDoc.shootProfile = integrationSweep;
        phDoc.shootDigitalBiomass = dbm(1);
        phDoc.rootDigitalBiomass = dbm(2);
        phDoc.rawSkeletonDistanceTransform = sampDist;
        phDoc.rawSkeletonDistanceTransformHistogram = SH;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate json string - OLD format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save json string - OLD format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make link to later
        JSON_string = savejson('carrotDoc',phDoc);
        % save JSON string
        filePushList{end+1} = [oPath n '_jdoc.j2son'];
        fileID = fopen(filePushList{end},'w');
        fprintf(fileID,strrep(JSON_string,'\/','\\/'));
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save json string - OLD format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % old JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % old csv file format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        % pet long path
        csvOutName = [oPath n '_petioleData.csv'];
        filePushList{end+1} = csvOutName;
        csvwrite(csvOutName,[petCount newPetCount petWidth petW mean(petL) petNUM petDIA]);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % old csv file format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push files to irods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(rPath)
            pushToiRods(rPath,filePushList);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push files to irods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      
    catch ME
        close all
        getReport(ME)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END and CLEAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    fprintf(['*************************************************************************************************\n']);
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
    fprintf(versionString);
    fprintf(['*************************************************************************************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % END and CLEAN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compile_directory = '/mnt/scratch1/phytoM/flashProjects/carrotPipeline/wholeCarrots/tmpSubmitFiles/';
    mkdir(compile_directory)
    CMD = ['mcc -d ' compile_directory ' -a im2single.m -m -v -R -singleCompThread singleWholeCarrotAnalyze.m'];
    eval(CMD);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check single image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileName = '/iplant/home/turnersarahd/carrotData/wholeData/03_10_2016_diallel/7155_rep1_15.13.46.NEF';
    fileName = '/iplant/home/turnersarahd/carrotData/wholeData/09_13_2015_diallel/50092_rep2_18.26.12.NEF';
    fileName = '/iplant/home/turnersarahd/carrotData/wholeData/09_15_2015_diallel/50150_rep1_15.59.14.NEF';
    fileName = '/iplant/home/turnersarahd/carrotData/wholeData/09_22_2016_mapping/60973_rep20_17.26.33.NEF';
    fileName = '/iplant/home/turnersarahd/carrotData/wholeData/09_22_2016_mapping/60974_rep18_17.09.50.NEF';
    %fileName = '/iplant/home/turnersarahd/carrotData/wholeData/09_13_2015_diallel/50090_rep3_18.00.32.NEF';
    singleWholeCarrotAnalyze(fileName,3350,200,40,'./output/','');
%}