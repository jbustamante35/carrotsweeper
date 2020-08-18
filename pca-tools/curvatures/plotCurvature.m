function plotCurvature(skel, k, c, m, smft, excl, sav, fnm, dout, fidx, dstr)
%% plotCurvature: display curvatures on mask
%
% Usage:
%    plotCurvature(skel, k, c, m, smft, excl, sav, fnm, dout, fidx, dstr)
%
% Input:
%   skel: post-processed binary mask
%   k: curvatures of each region
%   c: contours of each region
%   m: indices of maximum curvature for each region
%   smft: smoothing filter used for analysis
%   excl: number of columns excluded from the mask
%   fnm: filename for saving image
%   sav: boolean to save figure as png image
%   dout: path to save figure
%   fidx: index handle for figure
%
% Author Julian Bustamante <jbustamante@wisc.edu>
%

if nargin < 11
    dstr = 0; % Default to not showing distribution
    rows = 1;
    cols = 1;
end

%% Determine if data is split or unsplit
if isstruct(k.shoulder)
    splt  = 1;
    mclrs = {'r.' , 'm.' , 'c.'};
    mszs  = {1 , 15 , 15};
    ssz   = size(k.shoulder.upper, 1);
    tsz   = size(k.tip.upper, 1);
else
    splt = 0;
    ssz  = size(k.shoulder, 1);
    tsz  = size(k.tip, 1);
end

%%
figclr(fidx);
% lgnc = {'Contour' , 'Shoulders' , 'Tips' , 'Real Tip'};

% -------------------------------- Whole Mask -------------------------------- %
rgns  = {'whole' , 'shoulder' , 'tip'};
clrs  = {'b-' , 'g.' , 'y.'};
szs   = {3 , 8 , 8};

subplot(rows, cols, 1);
imagesc(skel);
hold on;

if splt
    % Plot contours of each region
    cellfun(@(x,y,z) plt(c.(x).upper, y, z), ...
        rgns, clrs, szs, 'UniformOutput', 0);
    cellfun(@(x,y,z) plt(c.(x).lower, y, z), ...
        rgns, clrs, szs, 'UniformOutput', 0);
    
    % Plot points of maximum curvature
    cellfun(@(x,y,z) plt(c.(x).upper(m.(x).upper,:), y, z), ...
        rgns, mclrs, mszs, 'UniformOutput', 0);
    cellfun(@(x,y,z) plt(c.(x).lower(m.(x).lower,:), y, z), ...
        rgns, mclrs, mszs, 'UniformOutput', 0);
else
    % Plot contours of each region
    cellfun(@(x,y,z) plt(c.(x), y, z), ...
        rgns, clrs, szs, 'UniformOutput', 0);
    
    % Plot points of maximum curvature
    cellfun(@(x,y,z) plt(c.(x)(m.(x),:), y, z), ...
        rgns, mclrs, mszs, 'UniformOutput', 0);
end

% legend(lgnc, 'Location', 'northeast');
ttl = sprintf('%s\n[%d columns removed | smooth filter %d]', ...
    fixtitle(fnm, 'carrots'), excl, smft);
title(ttl, 'FontSize', 10);

%% Plot Curvature Distributions
if dstr
    rows = 2;
    cols = 2;
    % --------------------------------- Shoulders -------------------------------- %
    subplot(rows, cols, 3);
    histogram(kS, 'EdgeColor', 'none', 'Normalization', 'pdf');
    ttl = sprintf('Shoulder Distributions [size %d]\nFilter size %d', ssz, smft);
    title(ttl, 'FontSize', 10);
    
    % ----------------------------------- Tips ----------------------------------- %
    subplot(rows, cols, 4);
    histogram(kT, 'EdgeColor', 'none', 'Normalization', 'pdf');
    ttl = sprintf('Tip Distributions [size %d]\nFilter size %d', tsz, smft);
    title(ttl, 'FontSize', 10);
end

% Save figure as png image
if sav
    saveFiguresJB(fidx, {fnm}, 0, 'png', dout);
end


end

