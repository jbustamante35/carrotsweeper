%%%%%%%%%%%%%%%%
% rotate trials to dims
function [t] = rotT(t)    
    % reform t
    d = 1:ndims(t);
    d = circshift(d,[0 1]);
    t = permute(t,d);
    t = r(t);
end
