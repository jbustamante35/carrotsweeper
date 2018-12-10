function ovec = getDateVec()
    dvec = datevec(now);
    ovec = [];
    for e = 1:numel(dvec)-1
        ovec = [ovec num2str(dvec(e)) '_'];
    end
    ovec(end) = [];
end