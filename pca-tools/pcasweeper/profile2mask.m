function msk = profile2mask(prf, wid)
%% profile2contour: plot a width profile as binary mask
% This function takes a width profile and generates a straightened mask
% representing that profile. The user can set the width of the mask to any
% number, but defaults to 300 pixels if left blank.
%
% Input:
%   prf: width profile
%   wid: arbitrary size to set mask width [default 300 pixels]
%
% Output:
%   msk: binary mask of the width profile
%

% Make blank mask template 2x the width and 1x the length of the profile
if nargin < 2
    % Need to choose an arbitrary width since these are normalized
    wid = 300; 
end

lng = size(prf,2);

% Convert profile to mask coordinates and cut in half
prf = prf * wid;
hlf = round(prf / 2);

% Generate mask by filling the template with ones the length of the profile
bot = zeros(size(hlf,1), lng);

for d = 1 : lng
    col = 1 : hlf(d);
    bot(col,d) = 1;
end

% Flip and Stitch
top = flipud(bot);
msk = [top ; bot];

end

