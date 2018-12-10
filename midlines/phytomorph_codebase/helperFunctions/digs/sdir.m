function [FileList] = sdir(cdir,type,Mpth)    
    % if iRODS via moniker then set string to '/'
    Mpth = attachSlash(Mpth);
    % make the list
    FileList = {};
    cnt = 1;
    for i = 1:size(cdir,1)
        if ~(cdir(i).isdir)
            [pth,nm,ext] = fileparts(cdir(i).name);
            if any(strcmp(type,ext(2:end)))
                FileList{cnt} = [Mpth cdir(i).name];
                cnt = cnt + 1;
            end
        end
    end
    FileList = FileList';
end