function [new_fileList] = convert2BSONList(fileList)
    import phytoCompute.*;
    % each set
    for s = 1:numel(fileList)
        % each file
        for e = 1:numel(fileList{s})
            % old file
            oldFile = fileList{s}{e};
            bFile = BSONfile(oldFile);
            new_fileList{s}{e} = bFile;
        end 
    end
end 