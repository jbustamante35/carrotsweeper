function [T] = generateDomainTransformation(PARA)
    LAMBDA = diag(PARA(2:3));
    T = [[cos(PARA(1)) -sin(PARA(1))];[sin(PARA(1)) cos(PARA(1))]];
    T = T*LAMBDA;
    %T = LAMBDA;
    T = [T ,[PARA(4);PARA(5)]];
end