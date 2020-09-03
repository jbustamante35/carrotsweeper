function [k , c , m , skel]  = computeCurvatures(msk, splt, cols2remove, flt, shoulder_size, tip_size)
%% computeCurvatures: compute curvatures of shoulders and tips
% This function takes a binary mask that can be straightened or unstraightened
% and computes the curvatures of the tip and the shoulder. The user can define
% the parameters for columns to exclude, smoothing curvature filter, and the
% number of coordinates to sample from the shoulder and tip.
%
% Explanation about splitting into upper and lower regions
%
% Usage:
%   [k , c , m , skel]  = ...
%        computeCurvatures(msk, splt, cols2remove, flt, shoulder_size, tip_size)
%
% Input:
%    msk: binary mask that can be straightened or unstraightened
%    splt: boolean to split into upper and lower regions
%    cols2remove: number of columns to exclude before sampling the shoulder
%    flt: smoothing filter size for the curvature function
%    shoulder_size: number of coordinates to sample from the shoulder
%    tip_size: number of coordinates to sample from the tip
%
% Output:
%    k: curvature values for each region
%    c: contour for each region
%    m: index of maximum curvature for each region
%    skel: processed mask after excluding columns
%
% Author Julian Bustamante <jbustamante@wisc.edu>
%

%% Default parameters
if nargin < 2
    splt                       = 1;
    [shoulder_size , tip_size] = deal(50);
    cols2remove                = 0;
    flt                        = 12;
end

%% Get curvatures and determine filter size and number of columns to remove
% Process Mask by removing left-most columns
[skel , cW] = maskProcessor(msk, cols2remove);

% Compute curvatures with range of smoothing filter size
mth      = 'closed';
[~ , kW] = cwtK(cW, flt, mth);
[~ , mW] = max(kW);

%% Split curvatures into shoulders and tips
% Shoulders are [(1 : n) , (N - n : n)]
cS = cW([1 : shoulder_size , end - shoulder_size + 1 : end], :);
kS = kW([1 : shoulder_size , end - shoulder_size + 1 : end], :);

% Tips are [(t - n) : t : (t + n)]
tCrd       = getTipIdx(skel);
[~ , tIdx] = min(cell2mat(arrayfun(@(x) pdist2(tCrd, cW(x,:)), ...
    1 : size(cW,1), 'UniformOutput', 0)'));

cT = cW(tIdx - tip_size : tIdx + shoulder_size, :, :);
kT = kW(tIdx - tip_size : tIdx + shoulder_size,:);

%% Split into top and bottom shoulder/tip or combine
flds = {'shoulder' ; 'tip' ; 'whole'};
if splt
    X       = {kS , kT , kW , cS , cT , cW};
    [Y , M] = cellfun(@(x) splitUpperAndLower(x), X, 'UniformOutput', 0);
    
    % Store regions into structure
    k = cell2struct(Y(1:3)', flds);
    c = cell2struct(Y(4:end)', flds);
    m = cell2struct(M(1:3)', flds);
    
else
    % Store regions into structure
    k = cell2struct({kS , kT , kW}', flds);
    c = cell2struct({cS , cT , cW}', flds);
    
    % Get max curvatures
    [~ , m.shoulder] = max(kS);
    [~ , m.tip]      = max(kT);
    m.whole          = mW;
end

end

function [skel , crv] = maskProcessor(msk, COLS2REMOVE)
%% maskProcessor: remove left-most columns from mask and get contour
% Description
%
% Usage:
%    [skel , crv] = maskProcesser(msk, COLS2REMOVE)
%
% Input:
%   msk: binary mask (right-left facing, black foreground, white background)
%   COLS2REMOVE: left-most columns to exclude from processing
%
% Output:
%   skel: processed mask with columns excluded
%   crv: contour corresponding to processed mask
%

%% Set constants for respective algorithms
if nargin < 2
    COLS2REMOVE = 15;  % Left-most columns to remove from the mask
end

% Other constants
MASK_THRESH     = 100;  % length to extend mask [removed after]
MIN_THRESH_SIZE = 100;  % number of columns to remove from contour [post-process]
FACE            = 3;    % Re-direct images left-right facing (original 3)
INTERP_SIZE     = 1000; % Number of coordinates to interpolate to

try
    %% Initial processing of the mask to face left-right and then pad
    if size(msk, 3) > 1
        msk = rgb2gray(msk);
    end
    
    % Force flip to left-right [arg = 3]
    %     skel = handleFLIP(msk, []);
    skel = handleFLIP(msk, FACE);
    chk  = ~imbinarize(skel);
    
    % Remove rows that are all empty
    while ~sum(chk(:,1))
        chk(:,1)    = [];
        COLS2REMOVE = COLS2REMOVE  + 1;
    end
    
    if COLS2REMOVE > 0
        skel = skel(:, COLS2REMOVE : end);
    end
    skel = ~padarray(skel, [0 MASK_THRESH], 'pre', 'replicate');
catch
    fprintf(2, '\nError with skeletonization\n');
    [skel , crv] = deal([]);
    return;
end

%% Post-processing of contour and mask
crv = getBWContour(skel, INTERP_SIZE);

% Remove padded area of mask and curve/midline coordinates
rmCrv        = crv(:,1) < MIN_THRESH_SIZE;
crv(rmCrv,:) = [];

% Remove padded area of mask
skel(:,1 : MIN_THRESH_SIZE) = [];
crv(:,1)                    = crv(:,1) - MIN_THRESH_SIZE;

% Re-interpolate contour
crv = interpolateOutline(crv, INTERP_SIZE);

end

function [y , m] = splitUpperAndLower(x)
%% splitUpperAndLower: split regions into upper and lower
sz = size(x);

switch mod(sz(1),2)
    case 0
        % Shoulders [half on each region]
        upper_range = 1                    : (sz(1) / 2);
        lower_range = upper_range(end) + 1 : sz(1);
    case 1
        % Tips [upper: half + tip | lower: 2nd half]
        upper_range = 1                    : (sz(1) / 2) + 1;
        lower_range = upper_range(end) + 1 : sz(1);
        
    otherwise
        % Something went wrong
        fprintf(2, 'Did something go wrong?\n');
        return;
end

switch sz(2)
    case 1
        % Curvatures: single vector
        y.upper = x(upper_range);
        y.lower = x(lower_range);
        
        % Get index of maximum curvature
        [~ , m.upper] = max(y.upper);
        [~ , m.lower] = max(y.lower);
    case 2
        % Coordinates: 2-dim vector
        y.upper = x(upper_range, :);
        y.lower = x(lower_range, :);
        m       = [];
    otherwise
        % Something went wrong
        fprintf(2, 'Did something go wrong?\n');
        y = [];
        return;
end

end
