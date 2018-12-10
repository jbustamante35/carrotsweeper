function [] = main_swell_HT(inFilePath,oPath,numberPerRow,numberToAnalyze,expectedImageNumber,comFunc,SKIP,template)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scan for new images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FileList = {};
    FileExt = {'tiff','TIF','tif'};
    verbose = 1;
    SET = sdig(inFilePath,FileList,FileExt,verbose);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% sort SET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        N = [];
        for img = 1:numel(SET{e})
            [p n ex] = fileparts(SET{e}{img});
            try
                N(img) = str2num(n);
            catch ME
                ME
            end
        end
        [N sidx] = sort(N);
        SET{e} = SET{e}(sidx);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove those which are in the junk folder
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        fidx = strfind(SET{e}{1},'junk');
        if ~isempty(fidx)
            rmidx(e) = 1;
        else
            rmidx(e) = 0;
        end
    end
    SET(find(rmidx)) = [];
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% find sets with more than Expected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        numImages(e) = numel(SET{e});
    end
    rmidx = comFunc(numImages,expectedImageNumber);
    SET(rmidx) = [];
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% generate output file names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        [swell{e} area{e} para{e} err{e} fit{e}] = generateOutFileBase(SET{e}{1},oPath,1);
    end
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove sets which have data in the oPath
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% remove compiled results
    oFileList = {};
    FileExt = {'csv'};
    verbose = 1;
    oSET = gdig(oPath,FileList,FileExt,verbose);
    rmidx = [];
    for e = 1:numel(oSET)
        fidx = strfind(oSET{e},'compiled_results');
        if ~isempty(fidx)
            rmidx(e) = 1;
        else
            rmidx(e) = 0;
        end
    end
    oSET(find(rmidx)) = [];
    %%% look for data which already is run through algo
    rmidx = zeros(numel(SET),1);    
    for e = 1:numel(oSET)
        fidx1 = strfind(oSET{e},filesep);
        fidx2 = strfind(oSET{e},'--');
        snipFile = oSET{e}(fidx1(end)+1:fidx2(end)-1);
        for i = 1:numel(swell)
            fidx1 = strfind(swell{i},filesep);
            fidx2 = strfind(swell{i},'--');
            isnipFile = swell{i}(fidx1(end)+1:fidx2(end)-1);
            if strcmp(snipFile,isnipFile)
                rmidx(i) = 1;
            end
        end
    end
    SET(find(rmidx)) = [];
    swell(find(rmidx)) = [];
    area(find(rmidx)) = [];
    para(find(rmidx)) = [];
    fit(find(rmidx)) = [];
    err(find(rmidx)) = [];
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% analyze the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tm = clock;
    mainBOX = 10^3*[1.1165 0.5835 8.9880 9.8160]; % fixed crop box due to template
    mainBOX = 1.0e+04*[0.106851000000000   0.051151000000000   0.909598000000000   1.039198000000000]; % new crop box for second widiv run
    %mainBOX = 10^4*[0.102051000000000   0.029551000000000   0.911998000000000   1.017598000000000];
    for e = 1:numel(SET)
        
                  
        [S A rowCount dB C] = measureStack_finalDeploy(SET{e},[500 500 500 500],0);
        
        
        
        
        % measure the data
        %[swellValue areaValue numberPerRow dB C] = measureStack(SET{e},[500 500 500 500],mainBOX,numberToAnalyze,numberPerRow,SKIP,0);
        
        
        % read the template file
        [p,nm,ext] = fileparts(SET{e}{1});
        templateFileName = [p filesep 'TemplateFile.csv'];
        template = csvread(templateFileName);
        template = template(:);
        SET{e}
        size(swellValue)
        % reshape the data
        sz2 = numberPerRow;
        sz1 = size(swellValue,1)/sz2;
        bS = reshape(swellValue',[size(swellValue,2) sz1 sz2]);
        bS = permute(bS,[2 3 1]);
        bA = reshape(areaValue',[size(areaValue,2) sz1 sz2]);
        bA = permute(bA,[2 3 1]);
        bS = reshape(bS,[size(bS,1)*size(bS,2) size(bS,3)]);
        bA = reshape(bA,[size(bA,1)*size(bA,2) size(bA,3)]);
        
        % template reshape
        template = template(:);
        UQ = unique(template);
        
        % loop over each row
        x0 = {};
        f = {};
        fixed = {};
        random = {};
        stats = {};
        UerrorInFit = {};
        
        %for genoType = 1:size(bS,1)
        for u = 1:numel(UQ)
            genoType = UQ(u);
            % get a unique genotype
            fidx = find(template==UQ(u));
            % get a row of data
            tmpData = squeeze(bS(fidx,:));
            tmpArea = squeeze(bA(fidx,:));
            % write to disk
            %csvwrite(strrep(swell{e},'#ROWNUM#',num2str(genoType)),tmpData');
            %csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');
            csvwrite(strrep(swell{e},'#ROWNUM#',num2str(genoType)),tmpData');
            csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');
            % remove those for which the method return 0 or if noise via
            % difference
          
            gtmpData = gradient(tmpData');
            nanidx = (abs(gtmpData') > .05) | (tmpData < -.05) | (tmpArea==0);
            nanidx = imdilate(nanidx,ones([1 4]));
            
            gidx = ~all(nanidx,2);
            tmpArea(nanidx) = NaN;
            tmpData(nanidx) = NaN;
            
            try
                % fit row of data
                [fixed{genoType},random{genoType},stats{genoType}] = fitSwellCurvePerGroup(tmpData(gidx,:));
            catch ME
                % fill in error and do not report the results for this
                % genotype
                gidx(genoType) = logical(0);
                ME
            end
            
            
            try                
                fidxl = 1:sum(gidx);
                Fidx = zeros(size(gidx));
                Fidx(gidx) = fidxl;
                for tr = 1:size(tmpData,1)
                    if gidx(tr)
                        x0{genoType}(tr,:) = fixed{genoType} + random{genoType}(:,Fidx(tr));
                        tm = 0:(size(tmpData,2)-1);
                        [f{genoType}((tr),:)] = swellCurveFit(x0{genoType}((tr),:),tm);
                        UerrorInFit{genoType}(tr) = nanmean(tmpData((tr),:) - f{genoType}((tr),:));
                    else
                        x0{genoType}(tr,:) = zeros(size(fixed{genoType}));
                        f{genoType}(tr,:) = zeros(1,size(tmpArea,2));
                        UerrorInFit{genoType}(tr) = NaN;
                    end
                end
            catch ME
                ME
            end
            
            
            %{
            % fit the data on kernel at a time            
            for tr = 1:size(tmpData,1)
                toFit = tmpData(tr,:);
                
                [x0{genoType}(tr,:) er(tr)] = fminsearch(@(X)mySwellFit(toFit,X),[10^4 .01]);                
                prediction_f{genoType}(tr,:) = func(x0{genoType}(tr,1),x0{genoType}(tr,2),1:3*size(tmpData,2));
                f{genoType}(tr,:) = func(x0{genoType}(tr,1),x0{genoType}(tr,2),1:size(tmpData,2));
                UerrorInFit{genoType}(tr) = mean(toFit - f{genoType}(tr,:));
                %{
                plot(f{genoType}(tr,:),'r');
                hold on
                plot(tmpData(tr,:),'b');
                plot(toFit,'b');
                drawnow                
                %}
            end
            %}
            %{
            hold off
            waitforbuttonpress
            %}
            try
                csvwrite(strrep(para{e},'#ROWNUM#',num2str(genoType)),[x0{genoType} tmpArea(:,1)]);
                csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');
                csvwrite(strrep(err{e},'#ROWNUM#',num2str(genoType)),UerrorInFit{genoType});
                csvwrite(strrep(fit{e},'#ROWNUM#',num2str(genoType)),f{genoType}');
            catch ME
                ME
            end
        end
    end
    fprintf(['Average time per stack:' num2str(etime(clock,tm)/numel(SET))/60 '\n']);
    
    
    
    
end


%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % single run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inFilePath = '/mnt/snapper/kernelSwellingData/Scott/rawData/';
    inFilePath = '/mnt/snapper/kernelSwellingData/Jeff/rawData/';
    oPath = '/mnt/snapper/kernelSwellingData/Jeff/return/';
    main_swell_HT(inFilePath,oPath,10,240,[],[],1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compile_directory = '/mnt/scratch1/phytoM/flashProjects/swelling/tmpSubmitFiles/';
    CMD = ['mcc -d ' compile_directory ' -m -v -R -singleCompThread main_swell_HT.m'];
    eval(CMD);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate master file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mPath = '/mnt/snapper/kernelSwellingData/Scott/return/';
    [MS rL] = createMasterFile(mPath);
    fileName = 'WIDIV_MasterSheet.csv';
    oPath = '/mnt/snapper/kernelSwellingData/Scott/';
    outName = [oPath fileName];
    cell2csv(outName,MS);


%}



