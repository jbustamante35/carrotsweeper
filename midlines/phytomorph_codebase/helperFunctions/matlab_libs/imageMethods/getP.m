function [parent] = getP(obj)
    classList = superclasses(obj);
    parent = classList;
    for e = 1:numel(classList)
        P = superclasses(classList{e});
        parent = setdiff(parent,P);
    end    
    if ~isempty(parent)
        parent = [parent{1} getP(parent{1})];
    end
end