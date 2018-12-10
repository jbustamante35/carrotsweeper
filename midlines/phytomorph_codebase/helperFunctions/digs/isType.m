function [fidx] = isType(FileList,Type)
    for i = 1:size(FileList,1)
        [pth,nam,type] = fileparts(FileList{i});
        if ~isempty(type)
            if any(strcmp(type,Type))
                fidx(i) = 1;
            else
                fidx(i) = 0;
            end
        else
            fidx(i) = 0;
        end
    end
    fidx = logical(fidx);
end