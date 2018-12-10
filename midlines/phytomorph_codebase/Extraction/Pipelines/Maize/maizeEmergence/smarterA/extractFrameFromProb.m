function [stepSig,frame] = extractFrameFromProb(prob,thresh,smallT,cAmount)
    sig = prob > thresh;
    v = sig(end);
    sig(end) = 0;
    sig = [zeros(size(sig));sig;zeros(size(sig))];
    sig = imclearborder(sig);
    sig = sig(2,:);
    sig(end) = v;
    sig = bwareaopen(sig,smallT);
    sig = imclose(sig,strel('disk',cAmount));
    
    
    % close only if last N frames are signaling emerge
    R = regionprops(logical(sig),'PixelIdxList');
    nsig = zeros(size(sig));
    for r = 1:numel(R)
        if any(R(r).PixelIdxList == numel(nsig))
            nsig(R(r).PixelIdxList) = 1;
        end
    end
    sig = nsig;
    
    
    
    fidx = find(sig);
    
    
    
    
    if ~isempty(fidx)
        frame = fidx(1);
    else
        frame = 0;
    end
    stepSig = sig;
end