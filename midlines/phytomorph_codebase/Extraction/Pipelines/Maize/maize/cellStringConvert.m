function [out] = cellStringConvert(in)
    if isa(in,'cell')
        out = cell2str(in);
    elseif isa(in,'char')
        out = str2cell(in);
    end
end

function [s] = cell2str(c)
    s = ',';
    for e = 1:numel(c)
        s = [s c{e} ','];
    end
end

function [c] = str2cell(s)    
    fidx = strfind(s,',');
    for e = 1:numel(fidx)-1
        c{e} = s(fidx(e)+1:fidx(e+1)-1);
    end
end


