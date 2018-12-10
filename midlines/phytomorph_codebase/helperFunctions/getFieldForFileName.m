function [ret] = getFieldForFileName(fileName,table,field)
    if ~iscell(fileName)
        fileName = {fileName};
    end
    fidx = zeros(numel(fileName),1);
    for f = 1:numel(fileName)
        fidx = fidx | strcmp(table.fileName,fileName{f});
    end
    
    ret = table.(field)(fidx);
end 