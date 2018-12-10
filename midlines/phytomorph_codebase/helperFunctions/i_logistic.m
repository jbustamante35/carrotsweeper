function [f] = i_logistic(x,L,k,xo)
    f = L*log(1+exp(k*(x-xo)));
end