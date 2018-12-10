function [S A] = measureStack(Stack,BOX,mainBOX,numtoMeasure,numCOLS,SKIP,disp)
    
    % create figures for display
    if disp
        h1 = figure;
        h2 = figure;
    else
        h1 = [];
        h2 = [];
    end
    
    % get the centers of the kernels
    C = getCenters(Stack{1},mainBOX);
    C = realignCenters(Stack{1},C,BOX);  
    
    
    % call to manual crop if centers are not found
    if isempty(C) | mod(size(C,1),numCOLS) ~= 0
        I = imread(Stack{1});
        [I mainBOX] = imcrop(I);
        C = getCenters(Stack{1},mainBOX);
    end
    
    % measure all if numtoMeasure is less than 0
    if numtoMeasure < 0
        numtoMeasure = numel(Stack);
    end
    
    
    % loop over the centers
    parfor e = 1:size(C,1)
        try
           %{
           [ST] = preLoadKernelImages(Stack,C(e,:),BOX,1,.50);
           [M] = hmm_ver0(ST);
            %}
           
           
           
           
           %uST = mean(double(ST),4);
           %[boundary currentArea MASK] = getKernelArea(uST,5000);
           %{
           sST = std(double(ST),1,4);
           toR = 50;
           R = double(reshape(ST(:,:,:,1:toR),[size(ST,1)*size(ST,2) size(ST,3) toR]));
           uD = mean(R,3);
           sD = std(R,1,3);
           for fr = 1:size(ST,4)
            
           end
           %}
           cnt = 1;
           tS = [];
           boundary = {};
           F = {};
           for img = 1:SKIP:numtoMeasure
                I = getKernelImage(Stack{img},C(e,:),BOX);
                %{
                if cnt == 1
                    tmpDB = [];
                    S = [];
                    tmpF = [];
                else
                    tmpDB = boundary{cnt-1};
                    tmpF = F{cnt-1};
                end
                %}
                %[boundary{cnt} currentArea F{cnt} S] = getKernelArea_ver2(I,10^4,tmpDB,tmpF,S,cnt,50);
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
                fprintf(['Done with kernel:' num2str(e) ':image:' num2str(img) '\n']);
           end
           
           A{e} = tS;
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

%{
    
%}