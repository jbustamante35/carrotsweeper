function [bi] = obtainBasisInformation_ver0(data)
    % decompose
    [S C U E L ERR LAM] = PCA_FIT_FULL_T(data,size(data,1));
    % store and return
    bi.E = E;
    bi.U = U;
end