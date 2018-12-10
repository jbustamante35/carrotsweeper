function [fidx] = watershedMax(in,para)
    % abs and flip for draining
    sig = -abs(in.K);
    % get watersheds
    L = watershed(sig);
    % find unqiue 
    UQ = unique(L);
    % get area and depth    
    for u = 1:numel(UQ)
        mask = L==UQ(u);
        sub_sig = sig(find(mask(:)));
        uH(u) = mean(sub_sig);
        uA(u) = sum(mask(:));
    end
    % calc min vol and find max volume of water shed
    VOL = sign(uH).*abs(uH).^.5.*uA;
    [uV sidx] = min(VOL);
    % find peak drain location
    mask = L==UQ(sidx);
    fidx = find(mask(1,:));
    values = sig(1,fidx);
    [peakV pidx] = min(values);
    fidx = fidx(pidx);
    
    
    %%% find point on the upper scale.
    %%% isolate region around it [2*PAD]
    %%% find max within pad
    %{
    WIDTH = 50;    
    % create global filter
    filter = ones(1,size(in.K,2));
    % loop over
    for i = 1:2
        % find the 
        midPoint = globalMNS(in.K(i,:).*filter,WIDTH);
        % create new filter
        F = zeros(size(filter));
        % create mask around max peak
        F(midPoint-WIDTH:midPoint+WIDTH) = 1;
        % refine filter
        filter = filter.*F;
    end
    % return
    fidx = midPoint;
    %}
end

function [midPoint] = globalMNS(sig,WIDTH)    
    % normalize the sig values
    nK = bindVec(-sig);
    % find peaks
    K = nonmaxsuppts(nK,40);    
    kidx = find(K);
    % find peaks heights
    samp = nK(kidx);
    % sample peaks heights
    [J selIDX] = max(samp);
    % return max peak
    midPoint = kidx(selIDX);
end