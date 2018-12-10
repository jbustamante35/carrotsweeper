function [fileList] = ffReplace(fileList,user,moniker)

    % 
    if nargin ~= 3
        % search for moniker for iRODS server
        moniker = '/phiRODS';
    end
    
    % 
    frontString = [filesep 'iplant' filesep 'home' filesep];
    
    % fuse file replace
    for ele = 1:numel(fileList)
        
        %
        oldFile = fileList{ele};
        % 
        fidx1 = strfind(oldFile,'/home/');
        fidx2 = strfind(oldFile,'/');
        if nargin == 1
            user = oldFile(fidx1+6:fidx2(3)-1);
        end
        % create temp string
        tmp_frontString = [frontString user filesep];
        
        % find the moniker
        fidx3 = strfind(oldFile,moniker);
        sidx = find(fidx2 > fidx3(1));
        % get the back string
        tmp_backString = oldFile(fidx2(sidx)+1:end);
        
        
        % make string
        newFile = [tmp_frontString tmp_backString];
        
        % new file
        fileList{ele} = newFile;
        
    end
end