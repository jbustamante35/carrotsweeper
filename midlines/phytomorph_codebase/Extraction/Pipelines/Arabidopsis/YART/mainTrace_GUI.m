function [] = mainTrace_GUI()
    % 1) Do you want to make cropboxes?
    
    % 2) If so, do you want to make another?
    
    % 3) Dark roots on light background or light roots on dark background?
    
    % 4) Do you want the choice to name each dataset?
    EXTRA = [200 50];
    %FileList = uipickfiles('FilterSpec','/home/nate/*');
    FileList = uipickfiles();
    
    if ~isempty(FileList)
        subCropBoxes = questdlg('Do you want to make crop boxes?','Sub Images?','Yes','No','Yes');

        rootColor = questdlg('Light roots on dark background or Dark roots on light background?','Light or Dark roots?','Light Roots','Dark Roots','Light Roots');

        toLAT = questdlg('Do you want to trace lat roots?','Lats?','Yes','No','Yes');
        
        subCropBoxes = strcmp(subCropBoxes,'Yes');
        rootColor = strcmp(rootColor,'Light Roots');
        toLAT = strcmp(toLAT,'Yes');
        
        oPath = uigetdir('','Where do you want to save results?');
        oPath = [oPath filesep];
        
        % obtain clicks for images
        for f = 1:numel(FileList)
            [X{f} BOX{f}] = handleSingleImage(FileList{f},subCropBoxes);
            [tp,outName{f},te] = fileparts(FileList{f});
        end     
        
        % compute engine
        % for each file
        for f = 1:numel(FileList)
            % for each crop box
            for cb = 1:numel(X{f})
                cnt = 1;
                
                % for each pair of points clicked on
                for pt = 1:2:size(X{f}{cb},1)
                    %{
                    %always trace lats
                    latTrace = 0;
                    if toLAT                        
                        if mod(pt,4) == 3
                            latTrace = 1;                            
                        end
                    end
                    %}
                    
                    % trace path via engine
                    %tmp_pheno = myTraceLow_engine(FileList{f},X{f}{cb}(pt:pt+1,:),rootColor,toLAT*latTrace,1,EXTRA);
                    [tmp_pheno mainTrace{f}{cb}] = myTraceLow_engine(FileList{f},X{f}{cb}(pt:pt+1,:),rootColor,1,1,EXTRA);
                    
                    % measure phenotype
                    
                    %{
                    para.SNIP = 5;
                    pheno{f}{cb}(cnt) = measurePheno(tmp_pheno,para);
                    clear tmp
                    if ~isempty(pheno{f}{cb}(cnt).LAT)
                        for e = 1:numel(pheno{f}{cb}(cnt).LAT)
                            tmp(e) = measurePheno(pheno{f}{cb}(cnt).LAT(e),para);
                        end
                        pheno{f}{cb}(cnt).LAT = tmp;
                    end
                    %}
                    cnt = cnt + 1;
                end
            end
        end
        
        % view/report engine        
        close all
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
        for f = 1:numel(FileList)
            I = imread(FileList{f});            
            for cb = 1:numel(X{f})                                  
                
                %phenoReport(pheno{f}{cb},oPath,outName{f},cb)
                
                [H1 H2]  = mainTrace{f}{cb}.measurePhenotypes();
                
                
                image(I);
                hold on
                %phenoView(pheno{f}{cb},1);
                mainTrace{f}{cb}.plot();
                if ~isempty(BOX{f})
                    axis([BOX{f}{cb}(1) BOX{f}{cb}(1)+BOX{f}{cb}(3) BOX{f}{cb}(2) BOX{f}{cb}(2)+BOX{f}{cb}(4)])
                end
                set(gca,'Position',[0 0 1 1]);
                axis off;
                imageName = [oPath outName{f} '_cropBox_' num2str(cb) '.tif'];
                
                
                
                saveas(gca,imageName);
                drawnow
                pause(1)
                
                
            end
        end
        
        
        
    end
    
end



function [X BOX] = handleSingleImage(fileName,subCropBoxes)    
    X = {};
    BOX = {};
    %%%%%%%%%%%%%
    % read the image
    I = imread(fileName);
    flag = 1;
    while flag
        %%%%%%%%%%%%%
        % assign toProcess if/not crop boxes
        if subCropBoxes
            close all
            imshow(I);
            hold on
            for e = 1:numel(BOX)
                rectangle('Position',BOX{e},'EdgeColor','r')
            end
            [toProcess BOX{end+1}] = imcrop();
            close all
        else 
            toProcess = I;
            clear I;
        end
        %%%%%%%%%%%%%
        % convert to gray scale if color
        if size(toProcess,3) == 3
            toProcess = rgb2gray(toProcess);
        end
        %%%%%%%%%%%%%
        % get pixel co-ordinates
        [x1 x2 V] = impixel(toProcess);
        %%%%%%%%%%%%%
        % displace to global frame
        if subCropBoxes
            x1 = x1 + BOX{end}(1);
            x2 = x2 + BOX{end}(2);
        end
        %%%%%%%%%%%%%
        % assign to output
        X{end+1} = [x1 x2];        
        %%%%%%%%%%%%%
        % ask if there are more crop boxes
        if subCropBoxes
            another = questdlg('Do you want to make another crop box?','Sub Images?','Yes','No','Yes');
            if strcmp(another,'No')
                flag = 0;
            end
        else
            flag =0;
        end
    end
end