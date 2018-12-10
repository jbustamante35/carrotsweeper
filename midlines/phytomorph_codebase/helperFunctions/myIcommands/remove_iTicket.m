function [FileList] = remove_iTicket(FileList)
    for e = 1:numel(FileList)
        ridx = strfind(FileList{e},'#');
        FileList{e} = FileList{e}(1:(ridx(1)-1));
    end
end