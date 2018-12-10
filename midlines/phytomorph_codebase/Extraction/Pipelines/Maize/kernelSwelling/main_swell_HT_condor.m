function [] = main_swell_HT_condor(inFilePath,oPath,numberPerRow,expectedImageNumber,numberToAnalyze,comFunc)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% analyze the data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tm = clock;
    mainBOX = 10^3*[1.1165 0.5835 8.9880 9.8160]; % fixed crop box due to template
    for e = 1:numel(SET)
        % measure the data
        [swellValue areaValue] = measureStack(SET{e},[800 800],mainBOX,numberToAnalyze,numberPerRow,0);
        
        
        % reshape the data
        sz2 = numberPerRow;
        sz1 = size(swellValue,1)/sz2;
        bS = reshape(swellValue',[size(swellValue,2) sz1 sz2]);
        bS = permute(bS,[2 3 1]);
        bA = reshape(areaValue',[size(areaValue,2) sz1 sz2]);
        bA = permute(bA,[2 3 1]);
        % loop over each row
        clear x0;
        clear f
        for genoType = 1:size(bS,1)
            % get a row of data
            tmpData = squeeze(bS(genoType,:,:));
            tmpArea = squeeze(bA(genoType,:,:));
            % write to disk
            csvwrite(strrep(swell{e},'#ROWNUM#',num2str(genoType)),tmpData');
            csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');
            try
                % fit row of data
                [fixed{genoType},random{genoType},stats{genoType}] = fitSwellCurvePerGroup(tmpData);
            catch ME
                ME
            end
            for tr = 1:size(tmpData,1)
                x0{genoType}(tr,:) = fixed{genoType} + random{genoType}(:,tr);
                tm = 0:(size(tmpData,2)-1);
                [f{genoType}(tr,:)] = swellCurveFit(x0{genoType}(tr,:),tm);
                UerrorInFit{genoType}(tr) = mean(tmpData(tr,:) - f{genoType}(tr,:));
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
            csvwrite(strrep(para{e},'#ROWNUM#',num2str(genoType)),[x0{genoType} tmpArea(:,1)]);
            csvwrite(strrep(area{e},'#ROWNUM#',num2str(genoType)),tmpArea');
            csvwrite(strrep(err{e},'#ROWNUM#',num2str(genoType)),UerrorInFit{genoType});
            csvwrite(strrep(fit{e},'#ROWNUM#',num2str(genoType)),f{genoType}');
        end
    end
    fprintf(['Average time per stack:' num2str(etime(clock,tm)/numel(SET))/60 '\n']);
    
    
    
    
end


%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inFilePath = '/mnt/snapper/kernelSwellingData/Scott/rawData/';
    oPath = '/mnt/snapper/kernelSwellingData/Scott/return/';
    main_swell_HT(inFilePath,oPath,120);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compile_directory = '/mnt/scratch1/phytoM/flashProjects/swelling/tmpSubmitFiles/';
    CMD = ['mcc -d ' compile_directory ' -m -v -R -singleCompThread main_swell_HT_condor.m'];
    eval(CMD);
%}


