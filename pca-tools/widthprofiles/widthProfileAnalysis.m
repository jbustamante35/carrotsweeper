function PW = widthProfileAnalysis(F, FPATHS, rootDir, numC, savpca, savcsv, vis, pcdim, outlier_pct, rng, fidxs)
%% widthProfileAnalysis: mutliple PCA analyses of width profiles
% Get PCA on Widths, Length-Width Normalized Profile, Tips, Shoulders.
% A neat pipeline to do all this shit from raw straightened masks. Output is in
% the form of PCA objects to all of these, and other data too but I'm not sure
% what else to store that can't be extrapolated from the PCA data.
%
% Usage:
%   PW = widthProfileAnalysis(F, FPATHS, rootDir, numC, savpca, savcsv, vis, ...
%                               pcdim, outlier_pct, rng, fidxs)
%
% Input:
%   F: structure containing widths prepared from prepWidthAnalysis function
%   FPATHS: file paths to images
%   rootDir: root directory to save Output 
%   numC: number of PCs to reduce width profiles to
%   savpca: boolean to save pca output in .mat files
%   savcsv: boolean to save output in .csv and .xls files
%   vis: visualize results
%   pcdim: PC score to sort by [default 1]
%   outlier_pct: percentage of top and bottom outliers to omit [default 0.05]
%   rng: number of profiles to visualize [default 5]
%   fidxs: figure handle to the 2 figures to generate [default 1:2]
%
% Output:
%   PW: PCA on width profiles

%
%% Default values
if nargin < 6
    pcdim       = 1;    % PC Score to sort outliers
    outlier_pct = 0.05; % Remove top and bottom outliers
    rng         = 5;    % Number of profiles to visualize in montage
    fidxs       = 1:2;  % Default figure handles
end

%% Prep timer and constants
sprA = repmat('-', 1, 80);
sprB = repmat('=', 1, 80);

tAll = tic;
fprintf('\n%s\nRunning Width Profile Analysis [Save PCA = %s | Save CSV = %s | Visualize = %s | Outlier Pct = %s | PC = %s]\n%s\n', ...
    sprB, num2str(savpca), num2str(savcsv), num2str(vis), ...
    num2str(outlier_pct), num2str(pcdim), sprA);

%% PCA on Width, Tips, and Shoulders [removing outliers]
nm      = strsplit(F.Name, '_');
nmidx   = strfind(F.Name, '_');
widnm   = F.Name(nmidx(end)+1:end);
W       = F.Profiles;
nrmL    = F.LengthNorm;
nrmW    = F.WidthNorm;
ttlWids = size(W, 1);
rgn     = nm{end};

t = tic;
% Widths
fprintf('Determined Normalization Method: [%s | %s lengths | %s widths]\n', ...
    rgn, nrmL, nrmW);
fprintf('Performing PCA on %d widths, shoulders, and tips...', ttlWids);
pcnm       = sprintf('%slength_%swidth_%s', nrmL, nrmW, rgn);

%% TODO: Replace this with the separate function instead
[PW, remW] = pcaOmitOutliers(W, numC, pcnm, outlier_pct, savpca, pcdim);
PATHFW     = FPATHS(remW);

%% Refresh some data after removing outliers
bakWids = ttlWids;
ttlWids = length(remW);

fprintf('Removed %d outliers from each dataset...', bakWids - ttlWids);
fprintf('DONE! [%.02f sec]\n', toc(t));

%% PCA on length-width normalized profiles 
%% [TODO]
t = tic;
fprintf('Performing PCA on %d orthonormalized profiles...', ttlWids);
fprintf('DONE! [%.02f sec]\n', toc(t));

%% Get Curvatures
t = tic;
fprintf('Computing curvatures of %d width profiles...', ttlWids);



fprintf('DONE! [%.02f sec]\n', toc(t));

%% Save Results in a CSV and .mat file
if savcsv
    t = tic;
    fprintf('Saving output in .csv and .xls files...');
    datDir = 'Output';
    outDir = sprintf('%s/%s', rootDir, datDir);
    saveCSV(PATHFW, PW, pcnm, ttlWids, outDir, pcdim);
    
    fprintf('...DONE! [%.02f sec]\n', toc(t));
end

%% Visualize results
if vis
    t    = tic;
    figs = fidxs;
    
    fprintf('Visualizing ranges of PC scores and masks...');
    fnms{1} = showRangeProfiles(PW, rng, widnm, pcdim, nrmL, nrmW, fidxs(1));
    fnms{2} = showRangeMasks(PW, rng, widnm, PATHFW, pcdim, nrmL, nrmW, fidxs(2));
        
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
function [p, remIdx] = pcaOmitOutliers(W, pcs, pcnm, outlier_pct, savpca, pcdim)
%% pcaOmitOutliers:
%
%
% Usage:
%
%
% Input:
%   W: vectorized width profiles
%   pcs: 
%   pcnm:
%   outlier_pct: 
%
% Output:
%

%% Initial PCA, Remove Outliers, Re-Run PCA
p           = pcaAnalysis(W, pcs, savpca, pcnm, 0);
[W, remIdx] = removeOutliers(W, p.PCAScores, outlier_pct, pcdim);
p           = pcaAnalysis(W, pcs, savpca, pcnm, 0);

end

function [W , remIdx] = removeOutliers(W, scrs, outlier_pct, pcdim)
%% removeOutliers: returns width profiles after removing outliers
outIdx = isoutlier(scrs(:,pcdim), 'percentiles', ...
    100 * [outlier_pct , 1 - outlier_pct]);
remIdx = find(~outIdx);
W      = W(remIdx,:);

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
strc = cell2struct([vals , scrs], flds, 2);
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

%% Save data in xls and csv files
% Eigenvectors, Means, and Curvatures
% Whole profile, Tips, and Shoulders

estr = sprintf('%s_%sVectors_PC%d', tdate, dataname, pcdim);

enm  = struct('EigVecs', evecs, 'Means', emns');
etbl = struct2table(enm);

etnm = sprintf('%s/%s.csv', ddir, estr);
writetable(etbl, etnm, 'FileType', 'text');

etnm = sprintf('%s/%s', ddir, estr);
writetable(etbl, etnm, 'FileType', 'spreadsheet');
end

% function vals = getNameID(FNAMES, id)
% %% getNameID: extract values from an id in a filename
% expr = sprintf('%s_(?<id>.*?)}', id);
% val  = regexpi(FNAMES, expr, 'names');
% vals = cellfun(@(x) char(x.id), val, 'UniformOutput', 0);

% end

% function [nrmL , nrmW] = getNormalization(W)
% %% getNormalization: determine normalization method used
% ZEROTHRESH = 0.2; % fewer than 20% should be 0 in normalized lengths
% zeroL      = numel(find(W == 0)) / numel(W); % Check if lengths have appended zeros
% maxW       = max(W, [], 'all'); % Check if width normalized to 1
% 
% if zeroL > ZEROTHRESH   
%     % Lengths have appended zeros [non-normalized]    
%     nrmL = 'original';
% else
%     % Lengths are normalized
%     nrmL = 'normalized';
% 
% end
% 
% if maxW > 1
%     % Widths are non-normalized
%     nrmW = 'original';
% else
%     % Widths are normalized
%     nrmW = 'normalized';
% end
% 
% end

function fnm = showRangeProfiles(pd, nRng, dnm, pcdim, nrmL, nrmW, fIdx)
%% showRangeProfiles: show range of profiles after sorting by PC scores
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
figclr(fIdx);
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

fnm = sprintf('%s_%s_%slength_%swidth_PC%d_%02dProfiles', ...
    tdate, dnm, nrmL, nrmW, pcdim, tot);
end

function fnm = showRangeMasks(pd, nRng, dnm, fpaths, pcdim, nrmL, nrmW, fIdx)
%% showRangeMasks: show range of masks after sorting by PC scores

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
figstr = cellfun(@(i,g,s,r,t) strcat(i,g,s,r,t), idxstr, str{1}, sepstr, ...
    str{2}, tabstr, 'UniformOutput', 0);
catstr = fixtitle(cat(2, figstr{:}), 'carrots');

% Get current version to determine which visualization to use
visver = getVisVersion;

% Show masks or profiles from sorted range
figclr(fIdx);
for i = 1 : numel(carIdx)
    idx     = carIdx(i);
    scr     = scrs(idx,:);
    msks{i} = scores2mask(scr, evecs, mns);
    
    % Show mask in individual subplots [ R2019a or lower ]
    if strcmpi(visver, 'sep')
        msk = scores2mask(scr, evecs, mns);
        subplot(nummsks, 1, i);
        myimagesc(msk);
        ttl = sprintf('%s %04d\n%s', dnm, idx, figstr{i});
        title(ttl, 'FontSize', 10);
    end
end

% Show montage in single subplot [ R2019b or higher ]
if strcmpi(visver, 'mon')
    buf = 5;
    montage(msks, 'Size', [1 , nummsks], 'BorderSize', [buf , buf]);
    axis off;
    ttl = sprintf('PC Score Range [%d %s | PC %d]\n%s', ...
        tot, dnm, pcdim, catstr);
    title(ttl, 'FontSize', 8);
end

fnm = sprintf('%s_%s_%slength_%swidth_PC%d_%02dMasks', ...
    tdate, dnm, nrmL, nrmW, pcdim, tot);
end

function clrmap = generateColorArray(itrs)
%% generateColorMatrix: generate a cell array of colors of specific size
clrs   = {'k', 'b', 'r', 'g', 'c', 'm'};
nReps  = ceil(itrs / numel(clrs));
clrmap = repmat(clrs, 1, nReps);
clrmap = clrmap(1 : itrs);

end

function v = getVisVersion
%% Get current version to determine which visualization to use
your_version = eval('version(''-release'')');
youryear     = str2double(your_version(end-2:end-1));
yourver      = your_version(end);
baseyear     = 19;
basever      = 'b';

% I get that there was probably a better way to do this but I don't have the
% time for that.
if strcmpi(basever, yourver)
    % Base version [b]
    if youryear >= baseyear
        % R2019b [ lowest year to use montage ] or greater
        v = 'mon';
    else
        % Lower than base year
        v = 'sep';
    end
else
    % Pre-release version [a]
    if youryear >= baseyear
        % R2020a [ lowest year to use montage ] or greater
        v = 'mon';
    else
        % Lower than base year [ includes R2019a ]
        v = 'sep';
    end
end
end

