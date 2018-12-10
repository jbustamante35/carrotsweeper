function [e] = maizeKN(file)
    if iscell(file)
        for i = 1:numel(file)
            e(i) = cnt(char(file{i}.getFullFileName()));
        end
    end
    if ischar(file)
        e = cnt(file);
    end
end

function [e] = cnt(file)
    fidx = strfind(file,'/');
    d = file(fidx(end-1)+1:fidx(end)-1);
    fidx = strfind(file,'_');
    e = numel(fidx)-1;
end