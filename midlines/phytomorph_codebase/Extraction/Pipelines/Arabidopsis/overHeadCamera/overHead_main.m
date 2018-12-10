function [] = overHead_main(sFileList,oPath,map,GMModel)
    try
        
        mkdir(oPath)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Starting image label and crop box detection.\n']);tic
        % label and find crop boxes
        [typeTable,cropTableCheckerBoard,cropTableRedSquares] = parseAndLabelStack(sFileList,map,[],[],false);
        fprintf(['Ending image label and crop box detection@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the file names for the day images only
        fileNames = getDayfileNames(typeTable);
        [pth,nm,ext] = fileparts(fileNames{1});
        
        %{
        fidx = strfind(nm,'_');
        nm = nm(1:fidx(1)-1);
        if ~(isdeployed())
            fidx = strfind(pth,filesep);
            pth = pth(fidx(end-1)+1:fidx(end)-1);
        end
        %}
        
        % get the names of the checker board images
        CheckerBoardfileNames = getCheckBoardfileNames(typeTable);
        % get the crop boxes for the checker board images
        checkerBoardBoxes = getFieldForFileName(CheckerBoardfileNames,cropTableCheckerBoard,'CropBoxes');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Starting checker board analysis.\n']);tic
        % getting checker board correction
        [tform,resEstimate] = checkerBoardAnalysis(CheckerBoardfileNames{1},checkerBoardBoxes{1},false);
        fprintf(['Ending image label and crop box detection@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Starting image label and crop box detection.\n']);tic
        % label and find crop boxes
        [typeTable,cropTableCheckerBoard,cropTableRedSquares] = parseAndLabelStack(sFileList,map,tform,resEstimate,false);
        fprintf(['Ending image label and crop box detection@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get crop boxes for the first N image
        fprintf(['Starting crop box merge analysis.\n']);tic
        N = 10;
        cropBOX = {};
        rm = [];
        for f = 1:N
            [cropBOX{f}] = getPlantBoxesForFile(cropTableRedSquares,fileNames{f});
            if numel(cropBOX{f}) ~= 9
                rm(f) = true;
            else
                rm(f) = false;
            end
        end
        cropBOX(find(rm)) = [];
        [cropBOX] = mergeCropBoxes(cropBOX);
        fprintf(['Found( ' num2str(numel(cropBOX)) ' ) crop-boxes. \n'])
        fprintf(['Ending crop box merge analysis@' num2str(toc) '\n']);
        %renderBoxesOnImage(fileNames{6},cropBOX,'r',3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Starting mask analysis.\n']);tic
        outTable = table;
        Area = [];
        for f = 1:numel(fileNames)
            I = double(imread(fileNames{f}))/255;
            [Area(f,:),QR_text{f}] = getPlantCell_and_Mask(I,cropBOX,GMModel,map);
        end
        fprintf(['Ending mask analysis@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather possible fields
        fprintf(['Starting QR analysis.\n']);tic
        FLDS = {};
        for f = 1:numel(QR_text)
            for b = 1:numel(QR_text{f})
                [flds] = mj_getFields(QR_text{f}{b});
                FLDS = [FLDS flds];
                FLDS = unique(FLDS);
            end
        end
        if isempty(FLDS);FLDS = {'NA'};end
        % get values
        HEADER = FLDS';
        for f = 1:numel(QR_text)
            fprintf(['Obtaining qr data from frame:' num2str(f) '\n']);
            for b = 1:numel(QR_text{f})
                fprintf(['Obtaining qr data from box:' num2str(b) '\n']);
                for h = 1:size(HEADER,1)
                    try
                        fprintf(['Obtaining header:' num2str(h) '\n']);
                        [~,value] = mj_getFields(QR_text{f}{b},HEADER{h,1});
                        if f == 1 && ~isempty(value{1})
                            HEADER{h,b+1} = value{1};
                        elseif f == 1 && isempty(value{1})
                            HEADER{h,b+1} = 'NA';
                        elseif  f ~= 1 && ~isempty(value{1}) && strcmp(HEADER{h,b+1},'NA')
                            HEADER{h,b+1} = value{1};
                        end
                    catch ME
                        getReport(ME)
                    end
                end                
            end
        end
        fprintf(['Ending QR analysis@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Starting output to disk.\n']);tic
        oArea = {};
        for tm = 1:size(Area,1)
            for p = 1:numel(cropBOX)
                oArea{p}(tm) = Area(tm,p);
            end
        end
        outTable = table(fileNames,oArea{1}(:),oArea{2}(:),oArea{3}(:),oArea{4}(:),oArea{5}(:),oArea{6}(:),oArea{7}(:),oArea{8}(:),oArea{9}(:));
        outTable = table2cell(outTable);
        outTable = [HEADER;outTable];
        %writetable(outTable,[oPath pth '__' nm '.csv']);
        oFile = [oPath pth '__' nm '.csv'];
        fprintf(['Writing to file-->' oFile '\n'])
        cell2csv(oFile,outTable);
        fprintf(['Ending output to disk@' num2str(toc) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        getReport(ME)
    end
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load from CyVerse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %websave(['./overHeadfuncToCallAPP.mat'],'https://de.cyverse.org/dl/d/7507C2DA-9F2F-41AB-8408-888670320E71/overHeadFuncAPP2.mat');
    %load(['./overHeadfuncToCallAPP.mat'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run single
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FilePath = '/mnt/tetra/nate/overHeadCrash/Run 4_104/';
    FilePath = '/mnt/tetra/nate/overHeadCrash/Run 3_204/';
    FileList = {};
    FileExt = {'jpg'};
    tic
    FileList = gdig(FilePath,FileList,FileExt,1);
    toc
    


%}