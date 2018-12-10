function [D,BV] = myGS(D,v)   
    % project 
    P = D*v;
    P = P*v';
    D = D - P;
    [J,D,J,BV,J,J,J] = PCA_FIT_FULL(D,(size(D,2)-1));
    %D = PCA_BKPROJ(C,BV,U);
end