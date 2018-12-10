function [sweepD] = sweepPCA(C,E,U,SD,compToSweep,np)
    u = mean(C,1);
    for c = 1:numel(compToSweep)
        sweepVec = linspace(-SD(compToSweep(c)),SD(compToSweep(c)),np);
        for e = 1:numel(sweepVec)
            tmp = u;
            tmp(compToSweep(c)) = sweepVec(e);
            sweepD(c,e,:) = PCA_BKPROJ(tmp,E,U);
        end
    end
end