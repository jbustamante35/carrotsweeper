function [localList] = xfer_get(remoteList,localPath,disjoint_union,disp,ticketList)
    try
        %%%%%%%%%%%%%%%%%%
        % xfer from irods OR idfunction for local
        % uses icommands to tranfer to local file system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if disjoint union is null
        % default to false
        if nargin == 2
            disjoint_union = logical(0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if ticket is given - give location for ticket
        if nargin <= 4
            ticketCMD = '';
        else
            ticketCMD = ' -t "<T>" ';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create template commands
        % hardcoded force overwrite
        cmd_template = ['iget -f -N 2 -V' ticketCMD ' "<S>" "<D>"'];
        %cmd_template = ['iget -f -V' ticketCMD ' "<S>" "<D>"'];
        local_template = [localPath '<D>'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the target directory
        mkdir(['.' localPath]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % report with wait bar
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if disp
            h = waitbar(0,'Transferring data from iRODS to local machine');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:numel(remoteList)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if it is a iRODS file then transfer to localPath
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isIRODS(remoteList{i})
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xfer file
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [pth,nm,ext] = fileparts(remoteList{i});
                cmd = strrep(cmd_template,'<S>',remoteList{i});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if disjoin union flag - then tag with i-number
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disjoint_union
                    cmd = strrep(cmd,'<D>',[localPath num2str(i) '-' nm ext]);    
                else
                    cmd = strrep(cmd,'<D>',[localPath nm ext]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if ticketCMD isempty
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isempty(ticketCMD)
                    cmd = strrep(cmd,'<T>',ticketList{i});
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % perform transfer
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['System call:' cmd '\n'])
                [status, result] = system(cmd,'-echo');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % report status of operation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if status == 0
                    fprintf(['xFer:success: ' [nm ext] '\n']);
                else
                    fprintf(['xFer:fail: ' [nm ext] '\n']);
                    fprintf(['command@' cmd '\n']);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get local file name for return
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disjoint_union
                    localList{i} = strrep(local_template,'<D>',[num2str(i) '-' nm ext]);
                else
                    localList{i} = strrep(local_template,'<D>',[nm ext]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if it is a localfile then do nothing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                localList{i} = remoteList{i};
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if display 
            if disp
                h = waitbar(i/size(remoteList,1),h,'Transferring data from iRODS to local machine');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % close wait bar
        if disp;close(h); end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % return the local list
        localList = localList';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        errorReport = getReport(ME);
        errorReport
        fprintf(['error@xfer_get\n']);
    end
end