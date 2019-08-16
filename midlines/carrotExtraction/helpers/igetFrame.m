function BV = igetFrame(crv, rad)
%% igetFrame: calculates the tangent and normal bundles for the curve
%
% Usage:
%   BV = igetFrame(crv, rad)
%
% Input:
%   crv: curve for analysis
%   rad: radius for pca
%
% Output:
%   BV: a sequence of basis vectors for the tangent and normal space
%

%%
BV = [];

for i = 1 : size(crv, 1)
    tBV = igetFrame_atP(crv, i, rad);
    BV  = cat(3, BV, tBV);   
end

% Make the first dim the trial number
BV = permute(BV,[3 1 2]);

end