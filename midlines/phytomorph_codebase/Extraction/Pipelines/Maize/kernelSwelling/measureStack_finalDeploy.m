function [S A rowCount dB C] = measureStack_finalDeploy(Stack,BOX,oPath,disp)
    mkdir(oPath);
    fileList = {};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort the filenames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(Stack)
        [p,nm,~] = fileparts(Stack{e});
        p = strrep(p,filesep,'_');
        num(e) = str2double(nm);
    end
    customName = p;
    [~,sidx] = sort(num);
    Stack = Stack(sidx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort the filenames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create figures for display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if disp
        h1 = figure;
        h2 = figure;
    else
        h1 = [];
        h2 = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create figures for display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the centers of the kernels, realign, and get row coun
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get centers
    [C,imageSZ,RxC,imageSZ,BOX] = getCenters(Stack{1});
    % align centers
    C = realignCenters(Stack{1},C,BOX);
    % get row count
    rowCount = RxC(1);
    % OLD row count
    %[rowCount] = countRows(imageSZ(1:2),C,100);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % view the centers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = imread(Stack{1});
    imshow(I,[]);hold on;plot(C(:,2),C(:,1),'r*')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call to manual crop if centers are not found
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(C) | mod(size(C,1),numCOLS) ~= 0
        I = imread(Stack{1});
        [I mainBOX] = imcrop(I);
        C = getCenters(Stack{1},mainBOX);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % measure all if numtoMeasure is less than 0
    if numtoMeasure < 0
        numtoMeasure = numel(Stack);
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over the centers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % image SKIP for debug
    SKIP = 1;
    % to measure for debug
    numtoMeasure = numel(Stack);
    % for each centroid
    parfor e = 1:size(C,1)
        try
            % count variable
            cnt = 1;
            tS = [];
            boundary = {};
            F = {};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for each image
            for img = 1:SKIP:numtoMeasure
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % try to read the image
                % if fail replace with previous image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                try
                    % get the bounding box from the image
                    I = getKernelImage(Stack{img},C(e,:),BOX);
                catch
                    % get the bounding box from the image
                    I = getKernelImage(Stack{img-1},C(e,:),BOX);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get kernel area and bounding box
                [boundary{cnt},currentArea] = getKernelArea(I,10^3,[.2 .5]);
                tS(cnt) = currentArea;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if disp
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disp
                    figure(h1);
                    imshow(I,[])
                    hold on
                    plot(boundary{cnt}(:,2),boundary{cnt}(:,1),'r');
                    hold off
                    figure(h2);
                    plot(calcPercentSwelling(tS));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % increment the counting variable
                cnt = cnt + 1;
                fprintf(['Done with kernel:' num2str(e) ':' num2str(size(C,1)) ':image:' num2str(img) '\n']);
            end
            % store the area for the e-th 
            A{e} = tS;
            % store the boundary for the e-th img
            dB{e} = boundary;
            fprintf(['Done with kernel:' num2str(e) '\n']);
        catch ME
            fprintf(['Error on kernel:' num2str(e) '\n']);
            A{e} = zeros(1,numel(Stack));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert area to matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = cell2mat(A');
    % convert to percent swelling
    % use the first three frames as init
    for e = 1:size(A,1)
        S(e,:) = calcPercentSwelling(A(e,:),3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kernel = 1:numel(dB)
        for img = 1:numel(dB{kernel})
            if ~isempty(dB{kernel}{img})
                dB{kernel}{img} = interp1(1:size(dB{kernel}{img},1),dB{kernel}{img},linspace(1,size(dB{kernel}{img},1),500));
            else
                dB{kernel}{img} = NaN*zeros(500,2);
            end
        end
    end
    
    for kernel = 1:numel(dB)
        tmpDoc = [];
        tmpDoc = generatePhenotypeNode(tmpDoc,S(kernel,:),{'percent_area'},'percent_area');
        tmpDoc = generatePhenotypeNode(tmpDoc,A(kernel,:),{'area'},'area');
        tmpDoc = generatePhenotypeNode(tmpDoc,C(kernel,:),{'centerPoint'},'centerPoint');
        tmpDoc = generatePhenotypeNode(tmpDoc,[C(kernel,:) C(kernel,:) + 2*BOX(1,[1 3])],{'boundingBox'},'boundingBox');
        for time = 1:numel(dB{kernel})
            tmpDoc = generatePhenotypeNode(tmpDoc,dB{kernel}{time},{'dimension','along'},['boundary_time' num2str(time)] );
        end
        phDoc(kernel) = tmpDoc;
    end
    JSON_string = savejson('swellDoc',phDoc);
               
    
    fileList{end+1} = [oPath customName '_jdoc.json'];
    fileID = fopen(fileList{end},'w');
    fprintf(fileID,strrep(JSON_string,'\/','\\/'));
    fclose(fileID);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % show image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    showI = imread(Stack{1});
    h = image(showI);
    hold on
    for e = 1:size(C,1)
        plot(C(e,2),C(e,1),'r*');
        %text(C(e,2)-50,C(e,1)-50,num2str(e),'BackgroundColor','w');
    end
    axis equal
    axis off
    saveas(gca,[oPath filesep customName  '_imageCheck.jpg']);
    csvwrite([oPath filesep customName  '_area.csv'],A);
    csvwrite([oPath filesep customName  '_percentAreaIncrease.csv'],S);
    close all;
end

%{
    FilePath = '/mnt/tetra/nate/20171027_Scanner1/';
    FileList = {};
    FileExt = {'tiff'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    [S A rowCount dB C] = measureStack_finalDeploy(FileList,[150 150 150 150],'./output/',0);



%}