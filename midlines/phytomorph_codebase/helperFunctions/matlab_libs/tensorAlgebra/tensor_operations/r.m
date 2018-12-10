%%%%%%%%%%%%%%%%
% reshape
function [d s] = r(d)
    s = size(d);
    %d = reshape(d,[prod(s(1:end-1)) s(end)])';
    d = reshape(d,[s(1) prod(s(2:end))])';
end