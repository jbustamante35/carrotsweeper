function [pw, ps, pt, pr, wids, lens] = widthProfileAnalysis(rootDir, maskDir, savpca, savcsv, vis, pcdim, outlier_pct, rng)
%% widthProfileAnalysis: mutliple PCA analyses of width profiles
% Get PCA on Widths, Length-Width Normalized Profile, Tips, Shoulders.
% A neat pipeline to do all this shit from raw straightened masks. Output is in
% the form of PCA objects to all of these, and other data too but I'm not sure
% what else to store that can't be extrapolated from the PCA data.
%
% Usage:
%   [pw, ps, pt, pr, wids, lens] = widthProfileAnalysis(rootDir, maskDir, ...
%                                  savpca, savcsv, vis, pcdim, outlier_pct, rng)
%
% Input:
%   rootDir: path to root directory of image sub-directories
%   mskDir: name of the sub-directory containing straightened-masks
%   savpca: boolean to save pca output in .mat files
%   savcsv: boolean to save output in .csv and .xls files
%   vis: visualize results
%   pcdim: PC score to sort by
%   outlier_pct: percentage of top and bottom outliers to omit
%   rng: number of profiles to visualize
%
% Output:
%   pw: PCA on width profiles
%   ps: PCA on shoulders
%   pt: PCA on tips
%   pr: PCA on width-length normalized profiles
%   wids: widths to convert from normalized to raw value
%   lens: lengths to convert from normalized to raw value
%
%% Default values
if nargin < 6
    pcdim       = 1;    % PC Score to sort outliers
    outlier_pct = 0.05; % Remove top and bottom outliers
    rng         = 5;    % Number of profiles to visualize from sorted range
end

%% Prep timer and constants
sprA = repmat('-', 1, 80);
sprB = repmat('=', 1, 80);

tAll = tic;
fprintf('\n%s\nRunning Width Profile Analysis [Save PCA = %s | Save CSV = %s | Vis = %s | Out = %s | PC = %s]\n%s\n', ...
    sprB, num2str(savpca), num2str(savcsv), num2str(vis), ...
    num2str(outlier_pct), num2str(pcdim), sprA);

%% Get Raw Width Profiles from Straightened Masks
t = tic;
fprintf('Loading width profiles from %s/%s...', rootDir, maskDir);

[DSTS , ~, FPATHS] = loadWidthProfiles(rootDir, maskDir);
ttlWids            = numel(DSTS);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Prepare Normalized Width Profiles and Split Tips and Shoulders
t = tic;
fprintf('Preparing %d width profiles and splitting shoulders/tips...', ttlWids);
FULLTRP              = 1000;
[CLIPSIZE , SPLTTRP] = deal(150);

% Full Widths
[FW , wids , lens] = prepWidthProfiles(DSTS, 'normalize', FULLTRP); % OG

% Shoulders [needs to be flipped]
FLPS = cellfun(@fliplr, DSTS, 'UniformOutput', 0);
SHLD = cellfun(@(x) x(1 : CLIPSIZE), FLPS, 'UniformOutput', 0);
SH   = fliplr(prepWidthProfiles(SHLD, 'normalize', SPLTTRP));

% Tips
TIPS = cellfun(@(x) x(1 : CLIPSIZE), DSTS, 'UniformOutput', 0);
TP   = prepWidthProfiles(TIPS, 'normalize', SPLTTRP);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% PCA on Width, Tips, and Shoulders [removing outliers]
t = tic;
fprintf('Performing PCA on %d widths, shoulders, and tips...', ttlWids);

[pcfw , pcsh , pctp] = deal(5);

% Widths
wfnm       = sprintf('%dNormalizedWidths', ttlWids);
[pw, remW] = pcaOmitOutliers(FW, pcfw, wfnm, ttlWids, outlier_pct, savpca, pcdim);
PATHFW     = FPATHS(remW);

% Shoulders
shnm       = sprintf('%dNormalizedShoulders', ttlWids);
[ps, remS] = pcaOmitOutliers(SH, pcsh, shnm, ttlWids, outlier_pct, savpca, pcdim);
PATHSH     = FPATHS(remS);

% Tips
tpnm       = sprintf('%dNormalizedTips', ttlWids);
[pt, remT] = pcaOmitOutliers(TP, pctp, tpnm, ttlWids, outlier_pct, savpca, pcdim);
PATHTP     = FPATHS(remT);

% Refresh some data after removing outliers
bakWids = ttlWids;
ttlWids = length(remW);
wids    = wids(remW,:);
lens    = lens(remW,:);

fprintf('Removed %d outliers from each dataset...', bakWids - ttlWids);
fprintf('DONE! [%.02f sec]\n', toc(t));

%% Graham-Schmidt Analyses to perform orthonormalization of width profiles
t = tic;
fprintf('Orthonormalization on %d width profiles...', ttlWids);

% Covariance in data using Canonical Correlations Analysis
COVARMETHOD = 'cca';
pscrs       = pw.PCAScores;
wids_lens   = [wids , lens];
[bb, ~]     = pcaRegression(wids_lens, pscrs, COVARMETHOD);

% Use regressors to remove length and width information from dataset
H = grahamSchmidt(bb(:,1), bb(:,2), pscrs);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% PCA on length-width normalized profiles
t = tic;
fprintf('Performing PCA on %d orthonormalized profiles...', ttlWids);

pcrw       = 5;
wrnm       = sprintf('%dNormalizedLengthWidth', ttlWids);
[pr, remR] = pcaOmitOutliers(H, pcrw, wrnm, ttlWids, outlier_pct, savpca, pcdim);
PATHPR     = FPATHS(remR);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Save Results in a CSV and .mat file
nmsW = 'Widths';
nmsS = 'Shoulders';
nmsT = 'Tips';
nmsR = 'Normalized';

if savcsv
    t = tic;
    fprintf('Saving output in .csv and .xls files...');
    datDir = 'Output';
    outDir = sprintf('%s/%s', rootDir, datDir);
    saveCSV(PATHFW, pw, nmsW, ttlWids, outDir, pcdim);
    saveCSV(PATHSH, ps, nmsS, ttlWids, outDir, pcdim);
    saveCSV(PATHTP, pt, nmsT, ttlWids, outDir, pcdim);
    saveCSV(PATHPR, pr, nmsR, ttlWids, outDir, pcdim);
    
    fprintf('...DONE! [%.02f sec]\n', toc(t));
end

%% Visualize results
if vis
    t = tic;
    fprintf('Visualizing ranges of PC scores...');
    
    figs    = 1 : 3;
    fnms    = cell(1, numel(figs));
    fnms{1} = showScoreRange(pw, rng, nmsW, pcdim, 1);
    fnms{2} = showScoreRange(ps, rng, nmsS, pcdim, 2);
    fnms{3} = showScoreRange(pt, rng, nmsT, pcdim, 3);    
    
    if savcsv
        fprintf('Saving %d figures...', numel(figs));
        saveFiguresJB(figs, fnms, 0, 'png', outDir);
    end
    fprintf('...DONE! [%.02f sec]\n', toc(t));
    
    t = tic;
    fprintf('Showing ranges of Masks...');
    fnms{1} = showProfileRange(pw, rng, nmsW, PATHFW, pcdim, 1);
    fnms{2} = showProfileRange(ps, rng, nmsS, PATHSH, pcdim, 2);
    fnms{3} = showProfileRange(pt, rng, nmsT, PATHTP, pcdim, 3);    
    
    if savcsv
        fprintf('Saving %d figures...', numel(figs));
        saveFiguresJB(figs, fnms, 0, 'png', outDir);
    end
    fprintf('...DONE! [%.02f sec]\n', toc(t));
    
end

fprintf('%s\nFinished Width Profile Analysis [%.02f sec]\n%s\n', ...
    sprA, toc(tAll), sprB);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, remIdx] = pcaOmitOutliers(W, pcs, pcnm, ttlWids, outlier_pct, savpca, pcdim)
%% pcaOmitOutliers:
%
%
% Usage:
%
%   
% Input:
%
%
% Output:
%

%% Initial PCA, Remove Outliers, Re-Run PCA
p           = pcaAnalysis(W, pcs, savpca, pcnm, 0);
[W, remIdx] = removeOutliers(W, p.PCAScores, ttlWids, outlier_pct, pcdim);
p           = pcaAnalysis(W, pcs, savpca, pcnm, 0);

end

function [W , remIdx] = removeOutliers(W, scrs, ttlWids, outlier_pct, pcdim)
%% removeOutliers: returns width profiles after removing outliers
[~, srtIdx] = sortrows(scrs, pcdim);
outTotal    = ceil(outlier_pct * ttlWids);
remIdx      = sort(srtIdx(1 + outTotal : end - outTotal));
W           = W(remIdx,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveCSV(FPATHS, P, dataname, ttlWids, outDir, pcdim)
%% saveCSV:
%
% saveCSV(PATHFW, pw, 'Widths', ttlWids, rootDir)
%

%% Get file paths [with outliers removed]
[FDIRS , FNAMES, ~] = cellfun(@(x) fileparts(x), FPATHS, 'UniformOutput', 0);

% Extract filename information
ids  = {'UID' , 'Genotype'};
vals = cellfun(@(x) getNameID(FNAMES, x), ids, 'UniformOutput', 0);

% Extract PC scores, EigenVectors, and Means
scrs  = P.PCAScores;
evecs = P.EigVecs;
emns  = P.MeanVals;

%% Store data in structure and table
flds = {'UID' , 'Genotype' , sprintf('%sPCs', dataname)};
strc = cell2struct([vals , scrs] , flds, 2);
tbl  = struct2table(strc);

% Save xls and csv files
dout = sprintf('%s_CarrotPCA_%dCarrots_%dGenotypes_PC%d_%s', ...
    tdate, ttlWids, numel(FDIRS), pcdim, dataname);
ddir = sprintf('%s', outDir);

if ~exist(ddir, 'dir')
    mkdir(ddir);
end

%% Save PC Scores
tnm1 = sprintf('%s/%s.csv', ddir, dout);
writetable(tbl, tnm1, 'FileType', 'text');

tnm2 = sprintf('%s/%s', ddir, dout);
writetable(tbl, tnm2, 'FileType', 'spreadsheet');

%% Save Eigenvectors and Means in xls and csv files [full, tips, shoulders]
estr = sprintf('%s_%sVectors_PC%d', tdate, dataname, pcdim);

enm  = struct('EigVecs', evecs, 'Means', emns');
etbl = struct2table(enm);

etnm = sprintf('%s/%s.csv', ddir, estr);
writetable(etbl, etnm, 'FileType', 'text');

etnm = sprintf('%s/%s', ddir, estr);
writetable(etbl, etnm, 'FileType', 'spreadsheet');
end

function vals = getNameID(FNAMES, id)
%% getNameID: extract values from an id in a filename
%
expr = sprintf('%s_(?<id>.*?)}', id);
val  = regexpi(FNAMES, expr, 'names');
vals = cellfun(@(x) char(x.id), val, 'UniformOutput', 0);

end

%%
function fnm = showScoreRange(pd, nRng, dnm, pcdim, fIdx)
%% showScoreRange: show range of curves after sorting by PC scores

% Extract info from the dataset
scrs = pd.PCAScores;
din  = pd.InputData;
dsim = pd.SimData;

tot           = length(scrs);
[dsrt , didx] = sortrows(scrs, pcdim);

% Create range of indices and color array
itr    = ceil(size(dsrt,1) / nRng);
cnt    = 1;
carIdx = didx([1 : itr : size(dsrt, 1) , size(dsrt,1)])';
clrs   = generateColorArray(numel(carIdx));
lgn    = cell(1, size(carIdx,1));

% A little trick to make sure concat has same digits as total
ndigs = num2str(numel(num2str(tot)));
dstr  = sprintf('%%0%sd', sprintf('%s', ndigs));

% Plot input and simulated range
set(0, 'CurrentFigure', fIdx);
cla;clf;
for i = carIdx
    clr = clrs{cnt};
    idx = didx(i);
    plot(din(idx,:), 'LineStyle', '-', 'Color', clr, 'LineWidth', 2);
    hold on;
    g(cnt) = plot(dsim(idx,:), 'LineStyle', '--', 'Color', clr, 'LineWidth', 2);
    
    % A little trick to make sure concat has same digits as total
    cmd    = sprintf('sprintf(''%s'', idx)', dstr);
    idxstr = eval(cmd);
    
    lgn{cnt} = sprintf('%s %s', dnm, idxstr);
    cnt = cnt + 1;
end

% Remove simulated curves and show legend and title
arrayfun(@(x) set(get(get(x, 'Annotation'), 'LegendInformation'), ...
    'IconDisplayStyle', 'off'), g, 'UniformOutput', 0);
legend(cat(1, lgn{:}), 'Location', 'bestoutside');

ttl = sprintf('PC Score Range [%d %s | PC %d]', tot, dnm, pcdim);
title(ttl, 'FontSize', 10);

fnm = sprintf('%s_%s_sortedPCs_PC%d_%dCarrots', dnm, tdate, pcdim, tot);

end

function fnm = showProfileRange(pd, nRng, dnm, fpaths, pcdim, fIdx)
%% showProfileRange: show range of profilesafter sorting by PC scores

% Extract data from the range
scrs  = pd.PCAScores;
evecs = pd.EigVecs;
mns   = pd.MeanVals;
tot   = length(scrs);

% Create range of indices
[dsrt , didx] = sortrows(scrs, pcdim);
itr           = ceil(size(dsrt,1) / nRng);
carIdx        = didx([1 : itr : tot , tot]);
nummsks       = nRng + 1;
msks          = cell(1, nummsks);

% Get Genotype and Root Number
[~, fnames, ~] = cellfun(@(x) fileparts(x), fpaths(carIdx), 'UniformOutput', 0);
ids            = {'Genotype', 'Root'};
str            = cellfun(@(x) getNameID(fnames, x), ids, 'UniformOutput', 0);

% Fix Genotype-Root Number string
spc    = 5;
idxstr = arrayfun(@(i) sprintf('(%s)+', num2str(i)), carIdx, 'UniformOutput', 0);
sepstr = cellstr(repmat('_', numel(fnames), 1));
tabstr = cellstr(repmat(repmat('+', 1, spc), numel(fnames), 1));
figstr = cellfun(@(i,g,s,r,t) strcat(i,g,s,r,t), idxstr, str{1}, sepstr, str{2}, tabstr, 'UniformOutput', 0);
catstr = fixtitle(cat(2, figstr{:}), 'carrots');

% Show masks or profiles from sorted range
set(0, 'CurrentFigure', fIdx);
cla;clf;
for i = 1 : numel(carIdx)
    idx     = carIdx(i);
    scr     = scrs(idx,:);
    msks{i} = scores2mask(scr, evecs, mns);
    
    % Show mask in individual subplots
    %     msk = scores2mask(scr, evecs, mns);
    %     subplot(nRng, 1, i);
    %     myimagesc(msk);    
    %     ttl = sprintf('%s %04d\n%s', dnm, cIdx, figstr{i});
    %     title(ttl, 'FontSize', 10);
end

% Show montage in single subplot
buf = 5;
montage(msks, 'Size', [1 , nummsks], 'BorderSize', [buf , buf]);
axis off;
ttl = sprintf('PC Score Range [%d %s | PC %d]\n%s', tot, dnm, pcdim, catstr);
title(ttl, 'FontSize', 8);

fnm = sprintf('%s_%s_sortedProfiles_PC%d_%02dCarrots', dnm, tdate, pcdim, tot);
end

function clrmap = generateColorArray(itrs)
%% generateColorMatrix: generate a cell array of colors of specific size
clrs   = {'k', 'b', 'r', 'g', 'c', 'm'};
nReps  = ceil(itrs / numel(clrs));
clrmap = repmat(clrs, 1, nReps);
clrmap = clrmap(1 : itrs);

end

