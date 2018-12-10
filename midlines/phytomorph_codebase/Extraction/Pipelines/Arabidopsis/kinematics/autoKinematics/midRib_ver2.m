function [midlineM] = midRib_ver1(FT,MSK,F1,F2)
    % find midline through mask
    WID = find(any(MSK,1));
    midline = find(any(MSK,2));
    midline = midline(1):midline(end);
    midline = sub2ind([size(F1,1) size(F1,2)],midline',round(mean(WID))*ones(numel(midline),1));
    midlineM = [];
    midlineM(:,1) = F1(midline);
    midlineM(:,2) = F2(midline);
   
    
    
end
   