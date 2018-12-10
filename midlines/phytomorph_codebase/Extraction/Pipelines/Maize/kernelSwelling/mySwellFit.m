function [error] = mySwellFit(data,X,TAU)
    if nargin == 2
        TAU = 1:size(data,2);
    end
    sim = func(X(1),X(2),TAU);
    error = norm(sim - data);
end

