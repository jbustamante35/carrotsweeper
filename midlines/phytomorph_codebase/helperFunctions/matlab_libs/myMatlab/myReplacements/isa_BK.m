function [r] = isa(obj,type)
    r = builtin('isa',obj,type);
    %{
    if iscell(type)
        r = wildCard(type);
        for e = 1:numel(type)
            r(e+1) = builtin('isa',obj,type{e});    
        end
        r = any(r);
    else
        r = builtin('isa',obj,type);
    end
    %}
end

function [r] = wildCard(type)
    r = any(strcmp(type,'*'));
end