function msk = profile2mask(pro)
%% profile2mask: generate a binary mask from a width profile
%
%
% Usage:
%   msk = profile2mask(pro)
%
% Input:
%   pro: width profile array
%
% Output:
%   msk: binary mask with 1 values set to widths
%

% Set all width values at first value == 0 or < 0 to 0
zIdx = find(pro <= 0);
if ~isempty(zIdx)
    pro(zIdx(1):end) = 0;
end

hlf    = ceil(pro / 2);
tmpmsk = zeros(max(hlf) , length(pro));

for i = 1 : length(hlf)
    p           = 1 : hlf(i);
    tmpmsk(p,i) = 1;
    %     tmpmsk(p,i) = hlf(i);
end

msk = [flipud(tmpmsk) ; tmpmsk];

end

