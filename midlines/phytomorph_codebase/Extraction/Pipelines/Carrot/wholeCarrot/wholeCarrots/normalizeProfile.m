function [d wd] = normalizeProfile(d,np)
    fidx = find(d==0);
    if ~isempty(fidx)
        ep  = (fidx(1)-1);
    else
        ep = numel(d);
    end
    d = interp1(1:ep,d(1:ep),linspace(1,ep,np));
    wd = d*max(d)^-1;
end