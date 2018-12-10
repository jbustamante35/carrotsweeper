function [pth] = attachSlash(pth)
    % if the path does not end in the slash - then attach    
    if ~any((strcmp(pth(end),'\') | strcmp(pth(end),'/')))
        if isIRODS(pth)
            pth = [pth '/'];
        else
            pth = [pth filesep];
        end
    end
end