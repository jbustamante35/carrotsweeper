function M = PCA_BKPROJ(C,E,U)
%%%%%%%%%%%%%%%%
% INPUTS:   C       : = components ("unique" fingerprint)M       : = data matrix
%           E       : = basis vectors 
%           U       : = mean of data
%%%%%%%%%%%%%%%%
% OUTPUTS:  M       : = create simulated data - from model -> create "data"
%%%%%%%%%%%%%%%%
% coeffs
M = (E*C')';                                          % project to get the coeffs
for i = 1:size(M,1)
    M(i,:) = M(i,:) + U;                                % subtract the mean
end