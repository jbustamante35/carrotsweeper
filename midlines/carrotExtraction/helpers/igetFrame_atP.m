function BV = igetFrame_atP(crv, idx, rad)
%% igetFrame_atP: twists curve vector at specific index
% Assume that arc-length is along dim1
%
% Usage:
%   BV = igetFrame_atP(crv, idx, rad)
%
% Input:
%   crv: 2D curve with dim1 along curve and dim2 imageDim
%   idx: index of position to obtain curve@
%   rad: local fallof for PCA
%
% Output:
%   BV:
%

%%
P     = crv(idx, :);
delta = bsxfun(@minus, P, crv);
delta = sum(delta .* delta, 2) .^ 0.5;
delta = delta < rad;
% delta = find(delta < rad); % Old method

%
seg   = crv(delta, :);
dseg  = diff(seg,1,1);
dseg  = mean(dseg,1);

%%
[~, ~, ~, BV, ~, ~, ~] = PCA_FIT_FULL(seg, size(seg,2));
BV(:,2)                = twistVec(BV(:,1));

%
if dseg * BV(:,1) < 0
    BV(:,1) = -BV(:,1);
end

%%
BV(:,2) = twistVec(BV(:,1));
BV      = flipud(BV);
end