function [C] = applyScannerNetwork(networkObject,v)






    [C] = PCA_REPROJ_T(v,networkObject.E,networkObject.U);
    if isdeployed
        C = networkObject.func(C);
    else
        C = networkObject.func(C);
        %C = networkObject.net(C);
    end
end



