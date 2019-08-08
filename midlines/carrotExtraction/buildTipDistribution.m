function [DOUT, crvs, xi, yi, vi] = buildTipDistribution(CNTRS, SMOOTH, BIN, MAG, OUTLIERS, MODE, sav, vis)
%% buildTipDistribution: initialize model to optimize tip-finding algorithm
%
%
% Usage:
%   [crvs_init, xinit, yinit, vinit] = buildTipDistribution(CRVSET, SMOOTH)
%
% Input:
%   CNTRS: cell array of contours to train initial distribution
%   SMOOTH: initial smoothing value to use for computing curvature
%   BIN:
%   MAG:
%   OUTLIERS:
%   MODE:
%   sav:
%   vis:
%
% Output:
%   DOUT: structure containing additional outputs to save in a .mat file
%   crvs: sorted array of all curvatures from training set of contours
%   xi:
%   yi:
%   vi:

try
    %% Determine if input are contours or curvatures
    sz = size(CNTRS{1}, 2);
    switch sz
        case 1
            % Input is in the form of curvatures (update model)
            crvs = cat(1, CNTRS{:});
        case 2
            % Combine and sort curvature profiles (initialize model)
            CRVS = cellfun(@(x) cwtK(x, SMOOTH, 'closed'), ...
                CNTRS, 'UniformOutput', 0);
            CRVS = cellfun(@(x) x.K, CRVS, 'UniformOutput', 0);
            
            % Curvature distribution with outliers removed
            crvs_all = cat(1, CRVS{:});
            crvs_all = sort(crvs_all, 'descend');
            nout     = round(length(crvs_all) * OUTLIERS);
            idx      = nout : (length(crvs_all) - nout);
            crvs     = crvs_all(idx);
        otherwise
            % Incorrect input, error out
            fprintf('Incorrect should be entered as a cell array\n');
            crvs = [];
    end
    
    %% Set left-right limits to capture main distribution around mean
    LMAX       = MAG * std(crvs) + median(crvs);
    [~ , LMAX] = min(abs(LMAX - crvs));
    
    LMIN       = -MAG * std(crvs) + median(crvs);
    [~ , LMIN] = min(abs(LMIN- crvs));
    
    fprintf('ida: %d | idb: %d | mag: %.06f \n', LMIN , LMAX, MAG);
    
    %% Generate curvature probability distribution
    [yi, xi] = hist(crvs, linspace(crvs(LMIN), crvs(LMAX), BIN));
    yi(1)    = [];
    yi(end)  = [];
    xi(1)    = [];
    xi(end)  = [];
    
    % Sum of -logs of curvature probabilities
    yi = yi/ sum(yi);
    vi = -log(yi);
    
    %% Store output in a structure and save
    Ncar = numel(CNTRS);
    Ncrv = length(crvs);
    DOUT = v2struct(Ncar, Ncrv, crvs, xi, yi, vi, LMAX, LMIN);
    if sav
        nm = sprintf('%s_CurveDistribution_%s_%dObjects_%dCurves', ...
            tdate('s'), MODE, Ncar, Ncrv);
        save(nm, '-v7.3', 'DOUT');
    end
    
    if vis
        plt([xi ; vi]', 'k-', 1);
        ttl = sprintf('Curvature Distribution -log(P(k))\n%d Carrots %d Curves', ...
            Ncar, Ncrv);
        title(ttl);
    end
    
catch e
    % Incorrect input
    fprintf('Error building model\n%s\n', e.getReport);
    [DOUT, crvs, xi, yi, vi] = deal([]);
end

end