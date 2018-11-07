function [ret] = strcontains(strarray,searchStr)
    for e = 1:numel(strarray)
        a = [];
        for s = 1:numel(searchStr)
            for sub = 1:numel(searchStr{s})
                if sub == 1;
                    a(s) = ~isempty(strfind(strarray{e},searchStr{s}{sub}));
                else
                    a(s) = a(s) | ~isempty(strfind(strarray{e},searchStr{s}{sub}));
                end
                
            end
        end
        ret(e) = all(a);
    end
end