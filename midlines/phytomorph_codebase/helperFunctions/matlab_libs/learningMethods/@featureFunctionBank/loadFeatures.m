function [] = loadFeatures(fileList,keyList)
    for e = 1:numel(fileList)
        [fileKey] = getFileKey(fileList{e});
        for key = 1:numel(keyList)
            outFile = [oPort fileKey '--' keyList{eky} '.mat'];
        end
    end

end
function [fileKey] = getFileKey(rawData)
    inFile = rawData{e}.fileName;
    [pth fileKey ext] = fileparts(inFile);
end