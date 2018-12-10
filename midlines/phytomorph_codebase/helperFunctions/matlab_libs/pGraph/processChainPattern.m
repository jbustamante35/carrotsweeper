function [y] = processChainPattern(p,x)
    % call process chain pattern
    % defas: 
    % p.c(e).f(unc)         := eth - function call
    % p.c(e).p(ara)         := eth - parameters for function-e
    % p.i(g).f(unc)         := gth - integrator function
    % p.i(g).p(ara)         := gth - parameters for integrator function
    %
    X{1} = x;
    for f = 1:numel(p.c)
        X{f+1} = p.c(f).f(X{f},p.c(f).p);
    end
    for g = 1:numel(p.i)
        y{g} = p.i(g).f(X,p.i(g).p);
    end
end
