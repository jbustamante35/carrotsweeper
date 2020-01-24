function [pw, ps, pt, pr, wids, lens] =  widthProfileAnalysis(rootDir, maskDir, savpca, savcsv, vis)
%% widthProfileAnalysis: mutliple PCA analyses of width profiles
% Get PCA on Widths, Length-Width Normalized Profile, Tips, Shoulders.
% A neat pipeline to do all this shit from raw straightened masks. Output is in
% the form of PCA objects to all of these, and other data too but I'm not sure
% what else to store that can't be extrapolated from the PCA data.
%
% Usage:
%   [pw, ps, pt, pr, wids, lens] =  widthProfileAnalysis(rootDir, maskDir, sav)
%
% Input:
%   rootDir: path to root directory of image sub-directories
%   mskDir: name of the sub-directory containing straightened-masks
%   savpca: boolean to save pca output in .mat files
%   savcsv: boolean to save output in .csv and .xls files
%
% Output:
%   pw: PCA on width profiles
%   ps: PCA on shoulders
%   pt: PCA on tips
%   pr: PCA on width-length normalized profiles
%   wids: widths to convert from normalized to raw value
%   lens: lengths to convert from normalized to raw value
%
%% Prep timer and constants
sprA = repmat('-', 1, 80);
sprB = repmat('=', 1, 80);

tAll = tic;
fprintf('\n%s\nRunning Width Profile Analysis [Save PCA = %s | Save CSV = %s | Vis = %s]\n%s\n', ...
    sprB, num2str(savpca), num2str(savcsv), num2str(vis), sprA);

%% Get Raw Width Profiles from Straightened Masks
t = tic;
fprintf('Loading width profiles from %s/%s...', rootDir, maskDir);

[DSTS , ~, FPATHS] = loadWidthProfiles(rootDir, maskDir);
ttlWids            = numel(DSTS);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Prepare Normalized Width Profiles and Split Tips and Shoulders
t = tic;
fprintf('Preparing %d width profiles and splitting shoulders/tips...', ...
    ttlWids);
FULLTRP              = 1000;
[CLIPSIZE , SPLTTRP] = deal(150);

% Full Widths
[FW , wids , lens] = prepWidthProfiles(DSTS, 'normalize', FULLTRP);

% Shoulders
SHLD  = cellfun(@(x) x(1  : CLIPSIZE), DSTS, 'UniformOutput', 0);
SH    = fliplr(prepWidthProfiles(SHLD, 'normalize', SPLTTRP));

% Tips
TIPS = cellfun(@(x) x(end-(CLIPSIZE -1): end), DSTS, 'UniformOutput', 0);
TP   = fliplr(prepWidthProfiles(TIPS, 'normalize', SPLTTRP));

fprintf('DONE! [%.02f sec]\n', toc(t));

%% PCA on Width, Tips, and Shoulders
t = tic;
fprintf('Performing PCA on %d widths, shoulders, and tips...', ttlWids);

[pcfw , pcsh , pctp, pcrw] = deal(5);

% Full Widths
wfnm = sprintf('%dNormalizedWidths', numel(DSTS));
pw    = pcaAnalysis(FW, pcfw, savpca, wfnm, 0);

% Shoulders
shnm = sprintf('%dNormalizedShoulders', numel(SHLD));
ps   = pcaAnalysis(SH, pcsh, savpca, shnm, 0);

% Tips
tpnm = sprintf('%dNormalizedTips', numel(SHLD));
pt   = pcaAnalysis(TP, pctp, savpca, tpnm, 0);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Graham-Schmidt Analyses to perform orthonormalization of width profiles
t = tic;
fprintf('Orthonormalization on width profiles...');

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

wrnm = sprintf('%dNormalizedLengthWidth', numel(DSTS));
pr   = pcaAnalysis(H, pcrw, savpca, wrnm, 0);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Save Results in a CSV and .mat file
if savcsv
    t = tic;
    fprintf('Saving output in .csv and .xls files...');
    
    %
    [FDIRS , FNAMES , ~] = cellfun(@(x) fileparts(x), FPATHS, 'UniformOutput', 0);
    
    % Extract filename information
    ids  = {'UID' , 'Genotype'};
    vals = cellfun(@(x) getNameID(FNAMES, x), ids, 'UniformOutput', 0);
    
    % Extract PC scores
    allPCA = [pw , pt , ps , pr];
    scrs   = arrayfun(@(x) x.PCAScores, allPCA, 'UniformOutput', 0);
    
    % Store data in structure and table
    flds = {'UID' , 'Genotype' , 'ShoulderPCs' , 'TipPCs' , ...
        'WidthPCs' , 'NormalizedPCs'};
    strc = cell2struct([vals , scrs] , flds, 2);
    tbl  = struct2table(strc);
    
    % Save xls and csv files
    dout = sprintf('%s_CarrotPCA_%dCarrots_%dGenotypes', ...
        tdate, numel(DSTS), numel(FDIRS));
    dnm  = 'Output';
    ddir = sprintf('%s/%s', rootDir, dnm);
    mkdir(ddir);
    
    tnm1 = sprintf('%s/%s.csv', ddir, dout);
    writetable(tbl, tnm1, 'FileType', 'text');
    
    tnm2 = sprintf('%s/%s', ddir, dout);
    writetable(tbl, tnm2, 'FileType', 'spreadsheet');
    
    % Save Eigenvectors and Means in xls and csv files [full, tips, shoulders]    
    eouts = {'Width' , 'Tip' , 'Shoulder'};
    estr  = cellfun(@(x) sprintf('%s_%sVectors', tdate, x), ...
        eouts, 'UniformOutput', 0);
    evecs = {pw.EigVecs , pt.EigVecs , ps.EigVecs};
    emns  = {pw.MeanVals , pt.MeanVals , ps.MeanVals};
    
    for e = 1 : numel(estr)
        enm = struct('EigVecs', evecs{e}, 'Means', emns{e}');
        etbl = struct2table(enm);
        
        etnm = sprintf('%s/%s.csv', ddir, estr{e});
        writetable(etbl, etnm, 'FileType', 'text');
        
        etnm = sprintf('%s/%s', ddir, estr{e});
        writetable(etbl, etnm, 'FileType', 'spreadsheet');
    end
        
    fprintf('...DONE! [%.02f sec]\n', toc(t));
end

%% Visualize results
if vis
    t = tic;
    fprintf('Visualizing ranges of PC scores...');
    
    figs    = 1 : 4;
    fnms    = cell(1, numel(figs));
    rng     = 5;
    fnms{1} = showScoreRange(ps, rng, 'Shoulders', 1);
    fnms{2} = showScoreRange(pt, rng, 'Tips', 2);
    fnms{3} = showScoreRange(pw, rng, 'Widths', 3);
    fnms{4} = showScoreRange(pr, rng, 'Normalized', 4);
    
    if savcsv
        for f = 1 : numel(figs)
            fprintf('Saving figure %02d...', f);
            fnm = sprintf('%s/%s', rootDir, fnms{f});
            saveas(figs(f), fnm, 'png');
        end
    end
    
    fprintf('...DONE! [%.02f sec]\n', toc(t));
    
end

fprintf('%s\nFinished Width Profile Analysis [%.02f sec]\n%s\n', ...
    sprA, toc(tAll), sprB);

end

function vals = getNameID(FNAMES, id)
%% getNameID: extract values from an id in a filename
%
expr = sprintf('%s_(?<id>.*?)}', id);
val  = regexpi(FNAMES, expr, 'names');
vals = cellfun(@(x) char(x.id), val, 'UniformOutput', 0);

end

%%
function fnm = showScoreRange(pd, nRng, dnm, fIdx)
%% showScoreRange: show range of curves after sorting by PC scores
% Extract info from the dataset
scrs = pd.PCAScores;
din  = pd.InputData;
dsim = pd.SimData;

tot = length(scrs);
col = 1;
[dsrt , didx] = sortrows(scrs, col);

% Create range of indices and color array
itr    = ceil(size(dsrt,1) / nRng);
cnt    = 1;
carIdx = [1 : itr : size(dsrt, 1) , size(dsrt,1)];
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

ttl = sprintf('PC Score Range [%d %s]', tot, dnm);
title(ttl);

fnm = sprintf('%s_%s_sortedPCs_%dCarrots', dnm, tdate, tot);

end

function clrmap = generateColorArray(itrs)
%% generateColorMatrix: generate a cell array of colors of specific size
clrs   = {'k', 'b', 'r', 'g', 'c', 'm'};
nReps  = ceil(itrs / numel(clrs));
clrmap = repmat(clrs, 1, nReps);
clrmap = clrmap(1 : itrs);

end

