function [fileName] = disIRODS(fileName)


    if isIRODS(fileName)
        url = true;
    else
        url = false;
    end
    fileList = {fileName};
    tmpPath = [tempname filesep];
    mkdir(tmpPath);
    
    [fileName,iticket] = stripiTicket(fileName);
    
    if ~isempty(iticket)
        fileList = xfer_get(fileList,tmpPath,0,0,{[fileName iticket]});    
    else
        fileList = xfer_get(fileList,tmpPath,0,0);    
    end
    
    
    fileList = fileList{1};
    fileName = fileList;  
end