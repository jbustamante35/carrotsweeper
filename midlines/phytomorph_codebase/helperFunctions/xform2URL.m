function [fileList] = xform2URL(fileList)
    fileTemplate = 'http://davos.cyverse.org/irods-rest/rest/fileContents#fileName#';
    ticketTemplate = '?ticket=#ticket#';
    for e = 1:numel(fileList)
        tmpFile = fileList{e};
        tidx = strfind(tmpFile,'#');
        if ~isempty(tidx)
            ticket = tmpFile(tidx(1)+1:tidx(2)-1);
            tmpFile(tidx(1):end) = [];
        end
        fileList{e} = [strrep(fileTemplate,'#fileName#',tmpFile)];
        if ~isempty(tidx)
            fileList{e} = [fileList{e} strrep(ticketTemplate,'#ticket#',ticket)];
        end
    end
end
        