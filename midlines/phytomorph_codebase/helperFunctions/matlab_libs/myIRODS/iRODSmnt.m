function [mntPoint,login,password] = iRODSmnt(mntPoint)
    try
        worked = 0;
        while ~worked
            % get username and password
            [login password] = logindlg('Title','Login Title');    
            % init the iRODS connection
            worked = myiinit(login,password);
        end

        % make the mount point directory
        mkdir(mntPoint)
        % create irods mount
        cmd = ['irodsFs ' mntPoint ' -o max_readahead=0'];
        %cmd = ['/home/nate/irodsFs -ometadatacachetimeout 3600 -opreloadblocks 5 ' mntPoint filesep];
        % mount the drive
        [r,e] = system(cmd);
    catch
        
    end
end