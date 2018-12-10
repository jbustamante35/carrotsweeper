function [D,BV] = myGS(D,v,toSubDIM)   
    % project 
    for e = 1:4
        D = D - (D*v)*v';
    end
    [J,D,U,BV,J,J,J] = PCA_FIT_FULL(D,(size(D,2)-1));
    %[J,C,U,BV,J,J,J] = PCA_FIT_FULL(D,(size(D,2)-toSubDIM));
    %D = PCA_BKPROJ(C,BV,U);
    
end