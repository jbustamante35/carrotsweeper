function [] = main_swell_HT_ver2(inFilePath,oPath,expectedImageNumber,trimNumber)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scan for new images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FileList = {};
    FileExt = {'tiff','TIF','tif',''};
    verbose = 1;
    SET = sdig(inFilePath,FileList,FileExt,verbose);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% sort SET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        N = [];
        for img = 1:numel(SET{e})
            [p n ex] = fileparts(SET{e}{img});
            N(img) = str2num(n);
        end
        [N sidx] = sort(N);
        SET{e} = SET{e}(sidx);
    end
    if ~isempty(trimNumber)
        for e = 1:numel(SET)
            SET{e} = SET{e}(1:trimNumber);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% find sets with less then expectedImageNumber
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        numImages(e) = numel(SET{e});
    end
    rmidx = numImages ~= expectedImageNumber;
    SET(rmidx) = [];
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% generate output file names
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(SET)
        [swell{e} area{e} para{e} err{e} fit{e} mat{e}] = generateOutFileBase(SET{e}{1},oPath,1);
    end
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
    toRemove = zeros(numel(SET),1);    
    for e = 1:numel(oSET)
        fidx1 = strfind(oSET{e},filesep);
        fidx2 = strfind(oSET{e},'--');
        snipFile = oSET{e}(fidx1(end)+1:fidx2(end)-1);
        rmidx = zeros(numel(SET),1);    
        for i = 1:numel(swell)
            fidx1 = strfind(swell{i},filesep);
            fidx2 = strfind(swell{i},'--');
            isnipFile = swell{i}(fidx1(end)+1:fidx2(end)-1);            
            if strcmp(snipFile,isnipFile)
                rmidx(i) = 1;            
            else
                rmidx(i) = 0;
            end
            toRemove(i) = toRemove(i) || rmidx(i);
        end
    end
    
    
    
        
    SET(find(toRemove)) = [];
    swell(find(toRemove)) = [];
    area(find(toRemove)) = [];
    para(find(toRemove)) = [];
    err(find(toRemove)) = [];
    fit(find(toRemove)) = [];
    mat(find(toRemove)) = [];
    
    
    disp = 1;
    toSave = 1;
    %{
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% analyze the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    tm = clock;
    %mainBOX = 10^3*[0.0589 0.0512 0.9612 1.3176]; % fixed crop box due to template
    for e = 1:numel(SET)
        % measure the data
        [swellValue areaValue RxC dB C] = measureStack_ver2(SET{e},-1,1,0);
        % reshape the data    
        bS = reshape(swellValue',[size(swellValue,2) RxC(1) RxC(2)]);
        bS = permute(bS,[2 3 1]);
        bA = reshape(areaValue',[size(swellValue,2) RxC(1) RxC(2)]);
        bA = permute(bA,[2 3 1]);
        bC = reshape(C',[size(C,2) RxC(1) RxC(2)]);
        bC = permute(bC,[2 3 1]);
        bdB = squeeze(reshape(dB',[size(dB,1) RxC(1) RxC(2)]));
       
        
        % loop over each row
        clear x0;
        clear f
        for genoType = 1:size(bS,1)
            % get a row of data
            tmpData = squeeze(bS(genoType,:,:));
            tmpArea = squeeze(bA(genoType,:,:));
            tmpCenters = squeeze(bC(genoType,:,:));
            tmpbdB = squeeze(bdB(genoType,:));
            
            if toSave
                % write to disk
                csvwrite(strrep(swell{e},'#ROWNUM#',num2str(genoType)),tmpData');
                csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');            
            end
            % fit the data on kernel at a time            
            for tr = 1:size(tmpData,1)
                toFit = tmpData(tr,:);
                [x0{genoType}(tr,:) er(tr)] = fminsearch(@(X)mySwellFit(toFit,X),[10^4 .01]);                
                prediction_f{genoType}(tr,:) = func(x0{genoType}(tr,1),x0{genoType}(tr,2),1:3*size(tmpData,2));
                f{genoType}(tr,:) = func(x0{genoType}(tr,1),x0{genoType}(tr,2),1:size(tmpData,2));
                UerrorInFit{genoType}(tr) = mean(toFit - f{genoType}(tr,:));
                if disp
                    plot(f{genoType}(tr,:),'r');
                    hold on
                    plot(tmpData(tr,:),'b');
                    plot(toFit,'b');
                    %waitforbuttonpress
                    drawnow                
                end
                
                if all(toFit(10:end) > 0)
                    badFlag(genoType,tr) = 0;
                else
                    badFlag(genoType,tr) = 1;
                end
            end
            if toSave
                csvwrite(strrep(para{e},'#ROWNUM#',num2str(genoType)),x0{genoType});
                csvwrite(strrep(err{e},'#ROWNUM#',num2str(genoType)),UerrorInFit{genoType});
                csvwrite(strrep(fit{e},'#ROWNUM#',num2str(genoType)),f{genoType}');
                save(strrep(mat{e},'#ROWNUM#',num2str(genoType)),'tmpCenters','tmpbdB');
            end
        end
    end
    fprintf(['Average time per stack:' num2str(etime(clock,tm)/numel(SET)/60) '\n']);
    %{
        for g = 1:numel(x0)
            keep = find(~logical(badFlag(g,:)));
            tmpD = x0{g}(keep,:);
            uF(g) = mean(tmpD(:,2));
            sF(g) = std(tmpD(:,2),1,1)*numel(keep)^-.5;
        end
    %}
end


%{
    inFilePath = '/mnt/snapper/kernelSwellingData/Scott/rawData/';
    oPath = '/mnt/snapper/kernelSwellingData/Scott/return/';


    inFilePath = '/mnt/snapper/Manfred/rawData/';
    oPath = '/mnt/snapper/Manfred/return/'
    main_swell_HT_ver2(inFilePath,oPath,288,288);




    inFilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/kernel_swellData/';
    oPath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/kernel_swellData/return/'
    main_swell_HT_ver2(inFilePath,oPath,288,288);
%}


