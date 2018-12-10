function [bool] = isFUSE(P,moniker)
    
    if nargin == 1
        % search for moniker for iRODS server
        moniker = '/phiRODS';
    end
    
    
    if numel(P) > numel(moniker)
        bool = strfind(P,moniker);
    else
        bool = 0;
    end
    
end