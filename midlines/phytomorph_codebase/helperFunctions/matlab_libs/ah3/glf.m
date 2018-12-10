function [g] = glf(x,p)
    % p(1) := upper bound
    % p(2) := lower bound
    % p(3) := growthrate
    % p(4) := power
    % p(5) := coeff
    % p(6) := max
    %g = p(1) + (p(2)-p(1))*(1+p(5)*exp(-p(3)*(x-p(6)))).^(p(4)^-1);
    g = (1 + exp(-p(1)*(x-p(2)))).^-p(3);
end
%{
p = [1 0 3 3 1 1];

p = [10 0 100];
x = linspace(-1,1,1000);
g = glf(x,p);
hold on
plot(x,g)
hold all
%}
