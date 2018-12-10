function [Mdata sz] = loadFeatures(fileList,sindex)
    % load features from list of many files one feature
    Mdata = [];
    for e = 1:numel(fileList)
        load(fileList{e},'data');
        % if index is provided the load the sub index
        if nargin == 2
            index = loadIndex(fileList{e});
            [~,sidx,~] = intersect(index,sindex);
            data = data(:,sidx);
        end
        if nargout == 2
            sz{e} = loadSize(fileList{e});
        end
        Mdata = [Mdata;data];
    end
end

function [fileKey] = getFileKey(rawData,e)
    if isa(rawData,'imageFile')
        inFile = rawData{e}.fileName;
    else
        inFile = rawData{e};
    end
    [pth fileKey ext] = fileparts(inFile);
end

function [index] = loadIndex(index)
    load(index,'index')
end

function [sz] = loadSize(sz)
    load(sz,'sz');
end