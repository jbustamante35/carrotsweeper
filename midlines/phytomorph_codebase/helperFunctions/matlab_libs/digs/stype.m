function [idx] = stype(cdir,type)
    cnt = 1;
    idx = [];
    for i = 1:size(cdir,1)
        if ~(cdir(i).isdir)
            [pth,nm,ext] = fileparts(cdir(i).name);
            idx(cnt) = any(strcmp(type,ext(2:end)));
        else
            idx(cnt) = 0;
        end
        cnt = cnt + 1;
    end
    idx = find(idx);
end