function [flds,value] = mj_getFields(text,query_f)
    flds = {};
    value = {''};
    d1 = strfind(text,'{');
    d2 = strfind(text,'}');
    for e = 1:numel(d1)
        snip = text(d1(e)+1:d2(e)-1);
        dl = strfind(snip,'_');
        flds{e} = snip(1:dl(1)-1);
        value{e} = snip(dl(1)+1:end);
    end
    if nargin == 2
        if ~isempty(flds)
            value = value(strcmp(flds,query_f));
            flds = flds(strcmp(flds,query_f));
        else
            flds = {query_f};
            value = {''};
        end
    end
end