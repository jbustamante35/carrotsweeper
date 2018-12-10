function [L] = struct2list(S,field)
    L = {};
    for i = 1:numel(S)
        L{i} = S(i).(field);
    end
    L = L';
end