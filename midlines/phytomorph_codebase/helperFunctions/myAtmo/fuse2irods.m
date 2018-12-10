function [fileList] = fuse2irods(fileList,user,moniker)
    %%%%%%%%%%%%%%%%%
    % moniker 
    if nargin == 2
        % search for moniker for iRODS server
        moniker = '/phiRODS';
    end
    %%%%%%%%%%%%%%%%%
    % common front string
    frontString = [filesep 'iplant' filesep 'home' filesep];
    %%%%%%%%%%%%%%%%%
    % fuse file replace
    % for each file set
    for s = 1:numel(fileList)
        % for each file in set
        for e = 1:numel(fileList{s})

            % old file
            oldFile = fileList{s}{e};
            
            if isFUSE(oldFile,moniker)
            
                % create temp string
                tmp_frontString = [frontString user filesep];

                % find the file sep
                sep_fidx = strfind(oldFile,filesep);
                
                % find the moniker
                moniker_fidx = strfind(oldFile,moniker);
                sidx = find(sep_fidx > moniker_fidx(1));

                % get the back string
                tmp_backString = oldFile(sep_fidx(sidx)+1:end);

                % make string
                newFile = [tmp_frontString tmp_backString];
            else
                % pass through
                newFile = oldFile;    
            end
            
            % assign new file
            fileList{s}{e} = newFile;
            
            % verbose
            fprintf([strrep(oldFile,'\','\\') ' ---> ' strrep(newFile,'\','\\') '\n']);
        end
    end
end 