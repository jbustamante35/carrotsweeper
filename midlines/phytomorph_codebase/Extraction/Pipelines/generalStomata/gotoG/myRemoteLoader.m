function [data] = myRemoteLoader(fileName,varargin)
    if isIRODS(fileName)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate a temp local location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpPath = ['/mnt/tetra/nate/tmp/' tempname filesep];
        mkdir(tmpPath);
    end
      
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % get the first one for pre-allocate
    %%%%%%%%%%%%%%%%%%%%%%%%
    [rawName,ticket] = stripiTicket(fileName);
    if isIRODS(fileName)
        fileList = xfer_get({fileName},tmpPath,0,0);
    else
        fileList{1} = fileName;
    end
    
    data = load(fileList{1},varargin{:});
end