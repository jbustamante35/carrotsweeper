function [out , K] = cwtK(J, smth, flt)
%% cwtK: wrapper for continuous wavelet transform to compute curvature
% This function is essentially a wrapper for a customized cwt algorithm, which
% is the continous 1-dimensional wavelet transform.
%
% Usage:
%   [out , K] = cwtK(J, smth, flt)
%
% Input:
%   J: input curve's x-/y-coordinates
%   smth: smoothing parameter for cwt function
%   flt: run algorithm on 'open' or 'closed' contour
%
% Output:
%   out: output structure containing information about wavelet transform
%   K: the curvature array
%
% Additional notes from Nathan:
%   MUST REMOVE THE BASE
%   FIND AND SUPPRESS

%% If 'flt' not selected, default to 'open' algorithm
if nargin < 2; smth = 0;      end
if nargin < 3; flt  = 'open'; end

out = struct('K', [], 'baseSize', [], 'J', [], 'Filter', flt);
switch flt
    case 'open'
        try
            [out , K] = runOpenContour(J, smth);
        catch re
            fprintf(2, 'Error in cwtK (open) [%s]\n%s\n', flt, re.getReport);
        end

    case 'closed'
        try
            [out , K] = runClosedContour(J, smth);
        catch re
            fprintf(2, 'Error in cwtK (closed) [%s]\n%s\n', flt, re.getReport);
        end

    otherwise
        %% Default to run on open contour
        try
            [out , K] = runOpenContour(J, smth);
        catch re
            fprintf(2, 'Error in cwtK [%s]\n%s\n', flt, re.getReport);
        end
end
end

function [out , K] = runOpenContour(J, smth)
%% runOpenContour: subfunction to run CWT on open contour [default]
% This is the default behavior for this function

%% Calculate curvature on open contour
d1X1 = cwt(J(:,1), smth, 'gaus1');
d1X2 = cwt(J(:,2), smth, 'gaus1');
d2X1 = cwt(J(:,1), smth, 'gaus2');
d2X2 = cwt(J(:,2), smth, 'gaus2');
K    = (d1X1 .* d2X2 - d1X2 .* d2X1) .* (d1X1.^2 + d1X2.^2) .^ -3/2;

%% Smooth data
K(:, 1:3 * max(smth))           = 0;
K(:, end - 3 * max(smth) : end) = 0;

%% outPut
out.J        = J;
out.baseSize = sum(J(:, 2) == 1);
out.K        = K;
end

function [out , K] = runClosedContour(J, smth)
%% runClosedContour: subfunction to run CWT on closed contour
% Run curvelet wave transform on closed loop

%% Calculate curvature around closed contour
K = zeros(length(J), numel(smth));
for e = 1 : numel(smth)
    h0 = fspecial('gaussian', [1 , 5 * smth(e)], smth(e));
    h0 = h0 / sum(h0);

    % Calculate curvature
    tmp = imfilter(J, h0', 'circular');
    d1X = gradient(tmp')';
    d2X = gradient(gradient(tmp'))';

    %% Smooth range
    K(:,e) = (d1X(:,1) .* d2X(:,2) - d1X(:,2) .* d2X(:,1)) .* ...
        (d1X(:,1).^2 + d1X(:,2).^2).^-3/2;
end

%% Output curvature
out.J        = J;
out.baseSize = sum(J(:, 2) == 1);
out.K        = K;
end