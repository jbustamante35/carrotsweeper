function [S A rowCount dB C] = measureStack(Stack,BOX,mainBOX,numtoMeasure,numCOLS,SKIP,disp)
    
    % create figures for display
    if disp
        h1 = figure;
        h2 = figure;
    else
        h1 = [];
        h2 = [];
    end
    
    % get the centers of the kernels
    [C imageSZ] = getCenters(Stack{1},mainBOX);
    
    
    
    
    C = realignCenters(Stack{1},C,BOX);
    
    
    [rowCount] = countRows(imageSZ(1:2),C,100);
    
    
    %{
    % call to manual crop if centers are not found
    if isempty(C) | mod(size(C,1),numCOLS) ~= 0
        I = imread(Stack{1});
        [I mainBOX] = imcrop(I);
        C = getCenters(Stack{1},mainBOX);
    end
    %}
    
    % measure all if numtoMeasure is less than 0
    if numtoMeasure < 0
        numtoMeasure = numel(Stack);
    end
    
    
    % loop over the centers
    parfor e = 1:size(C,1)
        try
           cnt = 1;
           tS = [];
           boundary = {};
           F = {};
           for img = 1:SKIP:numtoMeasure
                I = getKernelImage(Stack{img},C(e,:),BOX);
                [boundary{cnt} currentArea] = getKernelArea(I,10^4);
                tS(cnt) = currentArea;
                if disp
                    figure(h1);
                    imshow(I,[])
                    hold on
                    plot(boundary{cnt}(:,2),boundary{cnt}(:,1),'r');
                    hold off
                    figure(h2);
                    plot(calcPercentSwelling(tS));
                end
                cnt = cnt +1;
                fprintf(['Done with kernel:' num2str(e) ':' num2str(size(C,1)) ':image:' num2str(img) '\n']);
           end
           
           A{e} = tS;
           dB{e} = boundary;
           fprintf(['Done with kernel:' num2str(e) '\n']);
        catch ME
            fprintf(['Error on kernel:' num2str(e) '\n']);
            A{e} = zeros(1,numel(Stack));
        end
    end
    A = cell2mat(A');
    for e = 1:size(A,1)
        S(e,:) = calcPercentSwelling(A(e,:),3);
    end
    close all;
end