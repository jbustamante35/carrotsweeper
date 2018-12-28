function out = cwtK(J, para)
%% cwtK: wrapper for continuous wavelet transform to compute curvature
% This function is essentially a wrapper for a customized cwt algorithm, which
% is the continous 1-dimensional wavelet transform. 
%
% Usage:
%   out = cwtK(J, para)
%
% Input:
%   J: input curve's x-/y-coordinates
%   para: smoothing value for cwt function
%
% Output:
%   out: output structure containing information about wavelet transform
%
% Additional notes from Nathan:
%   MUST REMOVE THE BASE
%   FIND AND SUPPRES

%%
out = struct('K', [], 'baseSize', [], 'J', []);

try
    %% Calculate curvature
    d1X1 = cwt(J(:,1), para{1}, 'gaus1');
    d1X2 = cwt(J(:,2), para{1}, 'gaus1');
    d2X1 = cwt(J(:,1), para{1}, 'gaus2');
    d2X2 = cwt(J(:,2), para{1}, 'gaus2');
    K    = (d1X1 .* d2X2 - d1X2 .* d2X1) .* (d1X1.^2 + d1X2.^2).^-3/2;
    
    %% Filter out poor data
    K(:, 1:3 * max(para{1}))           = 0;
    K(:, end - 3 * max(para{1}) : end) = 0;
    
    %% outPort
    out.K        = K;
    out.baseSize = sum(J(:, 2) == 1);
    out.J        = J;
catch ME
    fprintf(2, 'Error in cwtK\n%s\n', ME.getReport);
end
end
