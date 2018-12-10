function [pheno] = measurePheno(pheno,para)
    for e = 1:numel(pheno)
        % make linear regression for root
        lin = linspace(0,1,size(pheno(e).gamma,2));
        [p1 s1 u1] = polyfit(lin',pheno(e).gamma(1,:)',1);
        [p2 s2 u2] = polyfit(lin',pheno(e).gamma(2,:)',1);
        pheno(e).straightFit.p1 = p1;
        pheno(e).straightFit.u1 = u1;
        pheno(e).straightFit.s1 = s1;
        pheno(e).straightFit.p2 = p2;
        pheno(e).straightFit.u2 = u2;
        pheno(e).straightFit.s2 = s2;
        pheno(e).straightFitValue = [polyval(p1,[0 1]',s1,u1) polyval(p2,[0 1]',s2,u2)];
        % measure length
        dL = diff(pheno(e).gamma,1,2);
        dL = sum(dL.*dL,1).^.5;
        pheno(e).length = sum(dL);
        % tip angle
        str = max(size(pheno(e).gamma,2)-para.SNIP+1,1);
        [S C U E L ERR LAM] = PCA_FIT_FULL(pheno(e).gamma(:,str:end)',1);
        pD = diff(pheno(e).gamma(:,str:end),1,2);
        pD = mean(pD,2);
        if pD'*E < 0
            E = -E;
        end        
        pheno(e).E = E;
        pheno(e).tipAngle = acos(E(1)/norm(E))*180/pi;
        % measure the skewAngle
        tmp = (pheno(e).gamma(:,end) - pheno(e).gamma(:,1));
        pheno(e).skewAngle = acos(tmp(1)/norm(tmp))*180/pi;
    end
end