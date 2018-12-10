function [percentRed ampSig] = processMiniStack(miniStack,miniMask,fFrames,lFrames,GT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new step 1
    close all
    clear f
    for e = 1:size(miniStack,4)
        f(:,:,:,e) = imfilter(miniStack(:,:,:,e),fspecial('disk',5),'replicate');
        f(:,:,:,e) = rgb2hsv(f(:,:,:,e));
    end
    %{
    dL = diff(bsxfun(@times,f(:,:,1,:),miniMask),1,4);
    dL = sum(sum(dL.^2,1),2);
    plot(squeeze(dL));
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new step 2
    %close all
    %imshow(f(:,:,1,1),[])
    %%
    %close all
    miniMask = bwlarge(miniMask);
    MM = imdilate(imclearborder(miniMask),strel('disk',31,0));
    %imshow((f(:,:,1,1) < .05).*MM,[]);
    %figure;imshow((miniStack(:,:,1,1)),[]);
    %figure;imshow((miniStack(:,:,1,1).*(f(:,:,1,1) < .05)),[]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    [CM ampSig] = getChipMask_E_ver0(f);
    for e = 1:size(miniStack,4)
        tmp = miniStack(:,:,:,e)/255;
        tmp = bsxfun(@times,tmp,MM);
        imshow(imresize(tmp,2),[])
        imshow(imresize(cat(3,CM(:,:,e),CM(:,:,1),MM),2),[])
        %imshow(CM(:,:,e).*MM,[])
        title(num2str(e))
        %if e > 79
        %    waitforbuttonpress
    % end
        drawnow
    end
    
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% new step 3
    close all
    [CM ampSig] = getChipMask_E_ver0(f);
    percentRed = (bsxfun(@times,CM,MM));
    
    percentRed = sum(sum(percentRed,1),2);
    percentRed = squeeze(percentRed);
    percentRed = percentRed / sum(MM(:));
    
end
%{
    %dL = abs(dL);
    %dL = imfilter(dL,fspecial('average',[5 1]),'replicate');
    %ampSig = imfilter(ampSig',fspecial('average',[5 1]),'replicate');
    ampSig = ampSig';
    
    fMean = mean(dL(1:fFrames));
    lMean = mean(dL((end-lFrames):end));
    sdL = sort(dL,'descend');
    lMean = mean(sdL(1:fFrames));
    %dL = dL.*ampSig';
    rawSignal = dL;
    
    if abs((lMean - fMean)) > .03;

        
        [dL var1] = subBL((dL),fFrames);
        [ampSig  var2]= subBL((ampSig),fFrames);
        
        
        
        rawSignal = ((dL.*ampSig));
        rawSignal = rawSignal/std(rawSignal(1:fFrames));
        %{
        %ampSig = abs(ampSig);
        %dL = abs(dL);
        %ampSig = ampSig - min(ampSig);
        %dL = dL - min(dL);
        
        %{
        if mean(dL) < 0
            dL = -dL;
        end
        %}
        thresh = 5*(var1*var2);

        %dL = abs(dL);
         
        slope = gradient(rawSignal);
        slope = imfilter(slope,fspecial('average',[21 1]),'replicate');
        %}
        %{
        %plot(dL,'k')
        %hold all
        %dL = imfilter(dL,fspecial('average',[7 1]),'replicate');
        %dL = imfilter(dL,fspecial('average',[7 1]),'replicate');
        %plot(dL,'r');
        %plot(ones(size(dL))*mean(dL(1:100)));
        %plot(ones(size(dL))*mean(dL(1:100))+5*std(dL(1:100),1,1));
        
        %}
        
        if max(dL) > .05
            thresh = 20;%mean(rawSignal(1:fFrames))+20*std(rawSignal(1:fFrames),1,1);
        else
            thresh = inf;
        end
       
        
    else
        thresh = inf;
        
    end
    
    binarySig = rawSignal > thresh;
    fidx = find(binarySig);
    if ~isempty(fidx)
        binarySig(fidx(1):end) = 1;
    end
    %figure
    %imshow(miniStack(:,:,:,fidx(1))/255,[])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure
    for e = 1:size(miniStack,4)
        imshow(miniStack(:,:,:,e)/255,[])
        if  e >= fidx(1)
            waitforbuttonpress
        end
        title(num2str(e))
        drawnow
    end
    %}
end
%}