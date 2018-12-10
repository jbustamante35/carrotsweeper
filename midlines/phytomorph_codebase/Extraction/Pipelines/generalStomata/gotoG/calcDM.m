function [DIS] = calcDM(angleValues)
    NF = size(angleValues,2);
    cPhi = 0;
    DIS = zeros([NF NF size(angleValues,1) size(angleValues,3)]);


    for m = 1:size(angleValues,3)
        for rho = 1:size(angleValues,1)
            for k1 = 1:(NF-1)
                for k2 = k1:(NF-1)
                    DIS(k1+1,k2+1,rho,m) = minPhaseDifference(cPhi,k1,angleValues(rho,k1+1,m),k2,angleValues(rho,k2+1,m),false);
                end
            end
            tmp = triu(DIS(:,:,rho,m));
            DIS(:,:,rho,m) = tmp + -1*sign(tmp').*abs(tmp');
        end
    end


    DIS = permute(DIS,[3 1 2 4]);
end