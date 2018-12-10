function [fileList] = irods2fuse(fileList,moniker)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % moniker 
    if nargin == 1
        % search for moniker for iRODS server
        moniker = '/phiRODS';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ret opath] = system('echo $HOME');
    opath(end) = [];
    newBase = [opath moniker];    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
            if isIRODS(oldFile)
                % find the file sep
                fidx = strfind(oldFile,filesep);
                
                % create temp string
                tmp_frontString = oldFile(1:fidx(4));

                % get the back string
                tmp_backString = oldFile(fidx(4):end);

                % make string
                newFile = [newBase tmp_backString];
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