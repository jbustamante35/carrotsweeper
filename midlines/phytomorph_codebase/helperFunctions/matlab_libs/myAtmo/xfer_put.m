function [] = xfer_put(sourcePath,targetPath)
    %if is file vs directory
    [pth,nm,ext] = fileparts(sourcePath);
    if ~isempty(nm) % is file
        %
        [nPath,~,~] = fileparts(targetPath);
        % make the target directory
        cmd = ['imkdir -p "' nPath '"'];
        [status results] = system(cmd);
        % push data over
        cmd = ['iput -f "' sourcePath '" "' targetPath '"'];
        [status results] = system(cmd);
    else % is directory
        if ~(isempty(targetPath))
            if isIRODS(targetPath)
                %%%%%%%%%%%%%%%%%%
                % xfer to irods
                cmd = ['imkdir -p "' attachSlash(targetPath) '"'];        
                [status results] = system(cmd);
                % store old path    
                p = pwd;
                cmd = ['cd ' sourcePath];
                [status results] = system(cmd);    
                cmd = ['cd ' sourcePath ';iput -r -f * "'  attachSlash(targetPath) '"'];
                [status results] = system(cmd);
                % restore path
                cmd = ['cd ' p];
                [status,resul] = system(cmd);
            else
                mkdir(targetPath);
                copyfile(sourcePath,targetPath,'f');
            end
        end
    end
end