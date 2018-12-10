function [BV] = igetFrame_from_point(curve,idx,radius)
    % assume that arc-length is along dim1
    P = curve(idx,:);
    delta = bsxfun(@minus,P,curve);
    delta = sum(delta.*delta,2).^.5;
    delta = find(delta < radius);
    seg = curve(delta,:);
    dseg = diff(seg,1,1);
    dseg = mean(dseg,1);
    [sS sC sU BV sL sERR sLAM] = PCA_FIT_FULL(seg,size(seg,2));
    if dseg*BV(:,1) < 0
        BV(:,1) = -BV(:,1);
    end
    BV(:,2) = [BV(2,1);-BV(1,1)];
end