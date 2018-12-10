function [FileList] = id2name(S)    
    for i = 1:size(S,2)        
        cmd = ['iquest --no-page "SELECT DATA_NAME,COLL_NAME WHERE DATA_ID = ''' S{i} '''"'];
        [status query] = system(cmd);
        tmpS = iquery2struct(query);
        FileList{i} = [tmpS(1).COLL_NAME '/' tmpS(1).DATA_NAME];
    end
end