function [oS] = name2id(S)
    % about: this function will take a list of files and IF
    % the filelist is from iRODS then it will convert to the unique
    % id of the file ELSEIF the file list is from a local directory
    % then it will return the filelist
    
    % get the unique "collection" names
    for j = 1:size(S,1)
        [pth{j} tmp_nm{j} tmp_ext{j}] = fileparts(S{j});
        nm{j} = [tmp_nm{j} tmp_ext{j}];
    end
    % for each collection - make request in bulk
    UQ = unique(pth);
    cnt = 1;
    for u = 1:size(UQ,2);
        % search collection u
        if isIRODS(UQ{u})
            cmd = ['iquest --no-page "SELECT DATA_ID,DATA_NAME,COLL_NAME WHERE COLL_NAME = ''' UQ{u} '''"'];    
            [status query] = system(cmd);
            tmpS = iquery2struct(query);
            for i = 1:size(tmpS,2)
                lok = strcmp(nm,tmpS(i).DATA_NAME) & strcmp(pth,tmpS(i).COLL_NAME);
                if sum(lok) == 1
                    oS{cnt} = tmpS(find(lok)).DATA_ID;
                    cnt = cnt + 1;
                end
            end
        else
            fidx = find(strcmp(UQ{u},pth));
            for i = 1:numel(fidx)
                oS{cnt} = S{fidx(i)};
                cnt = cnt + 1;
            end
        end
    end
    oS = oS';
end