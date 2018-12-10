function [fileList] = stackFileList(fileLists)
    cnt = 1;
    fileList = [];
    for s = 1:numel(fileLists)
        for e = 1:numel(fileLists{s})
            fileList{cnt} = fileLists{s}{e};
            cnt = cnt + 1;
        end
    end
end