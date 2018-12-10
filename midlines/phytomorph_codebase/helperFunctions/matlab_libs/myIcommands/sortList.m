function [L] = sortList(L)
    
    for i = 1:numel(L)
        [pth{i} nm{i} ext{i}] = fileparts(L{i});
        num_nm(i) = str2num(nm{i});
    end
    [J sidx] = sort(num_nm);
    L = L(sidx);
end