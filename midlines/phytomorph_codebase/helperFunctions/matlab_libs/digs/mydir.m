function [dirList] = mydir(path)
    
    % set to pwd as default
    if nargin == 0
        path = pwd;
    end
    
    % if found then iRODS dir - else local system dir
    if isIRODS(path)
        % get the contents of the directory
        cmd = ['ils "' path '"'];
        [status,results] = system(cmd);
        % find the new lines
        fidx_nl = strfind(results,sprintf('\n'));
        % for each line
        for ln = 1:(size(fidx_nl,2)-1)
            % strip out line
            line = results(fidx_nl(ln)+1:fidx_nl(ln+1)-1);
            line = strtrim(line);
            % check for collection moniker
            colSig = strcmp(line(1:3),'C- ');
            %store result
            if colSig
                dirList(ln).name = line(4:end);
                dirList(ln).isdir = logical(1);
            else
                dirList(ln).name = [attachSlash(path) line(1:end)];
                dirList(ln).isdir = logical(0);
                % get properies
                %{
                cmd = ['iquest "SELECT DATA_SIZE WHERE DATA_NAME = ''' dirList(ln).name ''' AND COLL_NAME = ''' path ''''];
                [status,sz] = system(cmd);
                eidx = strfind(size,'=');
                nidx = strfind(size,sprintf('\n'));
                dirList(ln).bytes = str2num(sz(eidx+1:nidx-1));
                %}
            end
        end
        dirList = dirList';        
    else
        dirList = dir(path);
        dirList(1:2) = [];
        dirList = genFullPath(dirList,path);
    end
end

%{
    cdir = mydir(ipwd);
    
%}