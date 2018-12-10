function [transform] = toFromaffine(transform)
    if size(transform,1) < size(transform,2)
        transform = [transform;[zeros(1,(size(transform,2)-1)) 1]];
    else
        transform(end,:) = [];
    end
end