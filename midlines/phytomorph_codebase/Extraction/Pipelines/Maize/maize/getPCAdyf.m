function [SIM] = getPCAdyf(S,X,IDX)
    %{
    for i = 1:numel(IDX)
        C(:,:,i) = PCA_REPROJ([X(:,IDX) X(:,IDX(i))],S(1).E,S(1).U);
        
    end
    C(isinf(C)|isnan(C)) = 0;    
    %}
    
    
    C1 = PCA_REPROJ([X(:,1) X(:,end) X(:,1)],S(1).E,S(1).U);
    C1 = C1*diag(diag(S(1).LAM).^-.5);
    C2 = PCA_REPROJ([X(:,1) X(:,end) X(:,end)],S(end).E,S(end).U);
    C2 = C2*diag(diag(S(end).LAM).^-.5);
    C1(isinf(C1)|isnan(C1)) = 0;
    C2(isinf(C2)|isnan(C2)) = 0;
    CT = cat(3,C1,C2);
    %CT = mean(CT,3);
    dCT = diff(CT,1,3);
    dL = linspace(0,1,numel(S));
    
    
    for tm = 1:numel(S)
        tmp = CT(:,:,1) + dL(tm)*dCT;
        tmp = tmp*diag(diag(S(tm).LAM).^.5);
        tmp = PCA_BKPROJ(tmp,S(tm).E,S(tm).U);
        SIM(:,tm)=  tmp(:,end);
    end
end