function [C] = tenFunc(T,para)
    %%%%%%%%%%%%%%%%%%%%
    % para{1} := polyNum vec
    % para{2} := polyDegree vec
    % para{3} := domain cell vec
    %%%%%%%%%%%%%%%%%%%%
    func = spap2(para{1},para{2},para{3},T);    
    C = func.coefs;
end