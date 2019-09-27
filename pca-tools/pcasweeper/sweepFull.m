function S = sweepFull(pStruct, pcs, stps, idx)
%% sweepFull: run a full PCA sweep through a dataset
%
% 
% Usage:
%   S = sweepFull(pStruct, pcs, stps, idx)
%
% Input:
%   pStruct: output from my custom PCA structure
%   pcs: number of pcs to iterate through
%   stps: total number of steps up and down to sweep
%   idx: index in dataset to sweep; omit to use the mean representation
%
% Output:
%   S: structure containing data of up-sweep and down-sweep
%

%% Extract data from PCA structure
mns  = pStruct.MeanVals;
eigs = pStruct.EigVectors;
scrs = pStruct.PCAscores;

% Allocate output array and set up-down functions
S     = repmat(struct('up', [], 'mean', [], 'down', []), ...
    numel(pcs), numel(stps));
upFn  = @(x,y) x + y;
dwnFn = @(x,y) x - y;

% Determine sweep function
if nargin < 4
    swpFn = @(p,s) pcaSweep(mns, eigs, scrs, p, upFn, dwnFn, s, 0);
else
    swpFn = @(p,s) pcaSweep(mns, eigs, scrs, p, upFn, dwnFn, s, 0, idx);
end

%%
for pc = pcs
    for stp = stps
        fprintf('PC %d | Step %d\n', pc, stp);
          [~, S(pc,stp)] = swpFn(pc, stp);
    end
end

end

