function [oS] = name2struct(S)
    % get the unique "collection" names
    for j = 1:size(S,1)
        [pth{j} nm{j}] = fileparts(S{j});
    end
    % for each collection - make request in bulk
    UQ = unique(pth);
    cnt = 1;
    S = [];
    for u = 1:size(UQ,2);
        % search collection u
        cmd = ['iquest --no-page "SELECT DATA_ID,DATA_NAME WHERE COLL_NAME = ''' UQ{u} '''"'];    
        [status query] = system(cmd);
        tmpS = query2struct(query);
        for i = 1:size(tmpS,2)
           S(cnt) = tmpS(i); 
        end
    end
end