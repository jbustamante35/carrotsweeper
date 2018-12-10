function [path] = convertPath(path)
    if size(path,1) == 2
        dP = diff(path,1,2);
        dP = bsxfun(@plus,dP,[1 1]');
        path = [1 3]*dP;
        path(path>5) = path(path>5) - 1;
    end
end