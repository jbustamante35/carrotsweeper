function [dC] = differentialCore(X,Y,mA,mB,coreA,coreB,Yexplain)
    [A,B,r,U,V,stats] = canoncorr(X,Y);
    
    for e = 1:size(A,2)
        if sign(mA(:,e)'*A(:,e)) == -1
            A(:,e) = -A(:,e);
            B(:,e) = -B(:,e);
            U(:,e) = -U(:,e);
            V(:,e) = -V(:,e);
        end
    end
    
    
    CORE = corr(U,Yexplain);
    
    dC = abs(CORE - coreA);
    
    
    
    dC = .5*(dC.*abs(CORE).^-1 + dC.*abs(coreA).^-1);
    
    %{
    bar([CORE(1,:);coreA(1,:)]');
    waitforbuttonpress
    close all
    %}
end