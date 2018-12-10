function [SIM1 SIM2 SIM3 SIM4 SIM5 SIM6] = getSIM(init,final,U,E,L)
    
    %CON = PCA_REPROJ([init final],E{end},U{end});
    CONi = PCA_REPROJ([init init],E{1},U{1});
    CONi = CONi*diag(diag(L{1}).^-.5);
    %CONi(isinf(CONi,2) = 0;
    CONi(isnan(CONi)) = 0;
    
    CONf = PCA_REPROJ([init final],E{end},U{end});
    CONf = CONf*diag(diag(L{end}).^-.5);
    CONf(isnan(CONf)) = 0;
    
    for tm = 1:numel(U)
        tmp1 = CONi*diag(diag(L{tm}).^.5);
        tmp2 = CONf*diag(diag(L{tm}).^.5);
        guess1 = PCA_BKPROJ(tmp1,E{tm},U{tm});
        guess2 = PCA_BKPROJ(tmp2,E{tm},U{tm});
        SIM1(:,tm) = guess1(:,2);
        SIM2(:,tm) = guess2(:,2);
        SIM3(:,tm) = mean([guess1(:,2),guess2(:,2)],2);
        SIM4(:,tm) = guess2(:,1);
    end
    SIM5 = fliplr(SIM4);
    
    
    CONi = PCA_REPROJ([SIM4(:,1) SIM4(:,1)],E{1},U{1});    
    CONi = CONi*diag(diag(L{1}).^-.5);
    CONi(:,2) = 0;
    for tm = 1:numel(U)
        tmp1 = CONi*diag(diag(L{tm}).^.5);
        guess1 = PCA_BKPROJ(tmp1,E{tm},U{tm});
        SIM6(:,tm) = guess1(:,2);
    end
    
end