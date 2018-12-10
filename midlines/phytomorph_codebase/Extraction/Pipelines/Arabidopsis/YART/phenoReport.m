function [] = phenoReport(pheno,oPath,fileName,boxNumber)
    TA = [pheno.tipAngle];
    SK = [pheno.skewAngle];
    LN = [pheno.length];
    HEADER = {'Tip Angle','Skew Angle','Length'};
    
    % along the root axis
    for e = 1:numel(pheno)
        tK = cwtK(pheno(e).gamma',{5});
        for s = 1:size(tK.K,2)
            KS{s,e} = tK.K(s);
        end
    end
    
    % for the co-ordinates of root
    for e = 1:numel(pheno)
        str = (e-1)*2+1;
        for i = 1:size(pheno(e).gamma,2)
            XY{i,str} = pheno(e).gamma(1,i);
            XY{i,str+1} = pheno(e).gamma(2,i);
        end
    end
    
    % for the scaler phenotypes
    OS = {};    
    outM = [TA;SK;LN]';
    for e = 1:numel(HEADER)
        OS{1,e} = HEADER{e};
    end    
    for e = 1:size(outM,2)
        for p = 1:size(outM,1)
            OS{p+1,e} = outM(p,e);
        end
    end
    
    
    
    outName = [oPath fileName '_cropBox_' num2str(boxNumber) '_phenoOut.csv'];
    cell2csv(outName,OS,',');
    outName = [oPath fileName '_cropBox_' num2str(boxNumber) '_Curvature.csv'];
    cell2csv(outName,KS,',');
    outName = [oPath fileName '_cropBox_' num2str(boxNumber) '_Position.csv'];
    cell2csv(outName,XY,',');
end