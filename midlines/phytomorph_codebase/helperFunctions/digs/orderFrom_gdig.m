function [FileList] = orderFrom_gdig(tmpList,FileList)
    pth = {};
    %%% sep into sets
    for i = 1:numel(tmpList)
        [pth{i} nm{i}] = fileparts(tmpList{i});
    end
    
    if ~isempty(pth)
        %%% sep into sets
        [UQ ia ic] = unique(pth);
        for u = 1:size(UQ,2)
            FileList{end+1} = tmpList(ic==u);
        end
    else
        FileList = {};
    end
end