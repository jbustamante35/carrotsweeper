function [g] = swellCurveFit(alpha,t)
    g = alpha(1)*(1-exp(-alpha(2)*t));
    %g = alpha(1) + alpha(2)*t + alpha(3)*t.^-1;
    %e = g - data;
    %e = norm(e);
end