function [index] = selectIndex_ver0(data,selectFunction,transformation)
    % select data which meets the selection feature criteria after being
    % transformed by the tranformation function
    index = loadIndex(data);
    data = loadData(data);
    data = transformation(data);
    sidx = selectFunction(data);
    index = index(sidx);
end

function [data] = loadData(data)
    if ~isnumeric(data)
        load(data,'data');
    end
end

function [index] = loadIndex(index)
    if ~isnumeric(index)
        load(index,'index');
    end
end