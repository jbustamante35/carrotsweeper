function [bool] = isIRODS(path,moniker)
    
    if nargin == 1
        % search for moniker for iRODS server
        moniker = '/iplant';
    end
    
    
    if numel(path) > numel(moniker)
        bool = strcmp(path(1:numel(moniker)),moniker);
    else
        bool = 0;
    end
    
end