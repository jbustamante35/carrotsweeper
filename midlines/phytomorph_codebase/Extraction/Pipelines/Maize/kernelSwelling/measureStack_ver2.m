function [S A RxC dB nC] = measureStack_ver2(Stack,numtoMeasure,imageSKIP,disp)
    
    % create figures for display
    h1 = figure;
    h2 = figure;
    
    % get the centers of the kernel
    % added for auto detection
    [C uBOX RxC] = getCenters(Stack{1});
    uBOX = uBOX*2;
    BOX = round([uBOX(1) uBOX(1) uBOX(2) uBOX(2)]/2);
    nC = realignCenters(Stack{1},C,BOX);
    
    
    %{
        % show centers for debug
        I = imread(Stack{1});
        imshow(I,[]);
        hold on    
        plot(C(:,2),C(:,1),'k*')
        plot(nC(:,2),nC(:,1),'r*')
        for e = 1:size(nC,1)
            text(nC(e,2),nC(e,1),num2str(e));
        end
    %}
    
    
    %{
    % call to manual crop if centers are not found
    if isempty(C) | mod(size(C,1),numCOLS) ~= 0
        I = imread(Stack{1});
        [I mainBOX] = imcrop(I);
        C = getCenters_ver2(Stack{1},mainBOX);
    end
    %}
    
    if numtoMeasure < 0
        numtoMeasure = numel(Stack);
    end
    
    threshValues = [.25 .8];
    %C = fliplr(C);
    % loop over the centers
    parfor e = 1:size(nC,1)
        try
           tS = [];
           boundary = {};
           img = 1;
           for imgV = 1:imageSKIP:numtoMeasure
               
                I = getKernelImage(Stack{imgV},nC(e,:),BOX);
                
                
                if disp
                    [boundary{img} currentArea MASK] = getKernelArea(I,10^3,threshValues);
                else
                    [boundary{img},currentArea,MASK] = getKernelArea(I,10^3,threshValues);
                end
                
                
                tS(img) = currentArea;
                if disp
                    figure(h1);
                    imshow(I,[])
                    hold on
                    plot(boundary{img}(:,2),boundary{img}(:,1),'r');
                    plot(boundary{1}(:,2),boundary{1}(:,1),'g');
                    hold off
                    
                    figure(h2);
                    plot(calcPercentSwelling(tS,5));
                end
                fprintf(['Done with kernel:' num2str(e) ':image:' num2str(imgV) ':' num2str(numtoMeasure) '\n']);
                img = img + 1;
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
    S = [];
    for e = 1:size(A,1)
        S(e,:) = calcPercentSwelling(A(e,:),3);
    end
    
    %{
    % quick test
    S = S';
    S = reshape(S,[size(S,1) 16 12]);
    S = permute(S,[2 3 1]);
    for e = 1:size(S,1)
        tmp = squeeze(S(e,:,:));
        ridx = find((tmp(:,end) < .05) > 0);
        tmp(ridx,:) = [];
        U = squeeze(mean(tmp,1));
        ST = squeeze(std(tmp,1,1))*size(tmp,1)^-.5;
        errorbar(U,ST);
        hold all
        title(num2str(e))
        waitforbuttonpress
        LEG{e} = num2str(e);
    end
    legend(LEG)
    %}
    
    close all;
end

%{
  


%}

%{
    
%}