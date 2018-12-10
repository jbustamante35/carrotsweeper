function [e d] = mySigmoid_ver0(x,para,Y)
    e = [];
    d = para(1)./(1+exp(-para(2)*(x-para(3))));
    if nargin == 3
        e = norm(Y - d);
    end
end