function [swell_fileOut,area_fileOut,para_fileOut,err_fileOut] = generateOutFileBase(inFile1,oPath,numDeep)
    [pth,nm,ext]= fileparts(inFile1);
    fidx = strfind(pth,filesep);
    fileName = pth(fidx(end-numDeep)+1:end);
    fileName = strrep(fileName,filesep,'--');
    swell_fileOut = [oPath 'swell--' fileName '--' '#ROWNUM#' '.csv'];
    area_fileOut = [oPath 'area--' fileName '--' '#ROWNUM#' '.csv'];
    para_fileOut = [oPath 'para--' fileName '--' '#ROWNUM#' '.csv'];
    err_fileOut = [oPath 'error--' fileName '--' '#ROWNUM#' '.csv'];
end