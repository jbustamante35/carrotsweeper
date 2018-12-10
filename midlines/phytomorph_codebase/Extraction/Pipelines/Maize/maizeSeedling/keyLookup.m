function [v] = keyLookup(fileName,k)
    try
        fidx = strfind(fileName,k);
        gidx = strfind(fileName,'}');
        sidx = find(gidx > fidx(1));
        v = fileName((fidx(1)+numel(k)+1):(gidx(sidx(1))-1));
    catch ME
        v = '';
    end
end