function [pw, ps, pt, pr, wids, lens] =  widthProfileAnalysis(rootDir, maskDir, sav)
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
%   sav: boolean to save output in .mat files and .csv file
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
fprintf('\n%s\nRunning Width Profile Analysis [Save = %s]\n%s\n', ...
    sprB, num2str(sav), sprA);

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
pw    = pcaAnalysis(FW, pcfw, sav, wfnm, 0);

% Shoulders
shnm = sprintf('%dNormalizedShoulders', numel(SHLD));
ps   = pcaAnalysis(SH, pcsh, sav, shnm, 0);

% Tips
tpnm = sprintf('%dNormalizedTips', numel(SHLD));
pt   = pcaAnalysis(TP, pctp, sav, tpnm, 0);

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
pr   = pcaAnalysis(H, pcrw, sav, wrnm, 0);

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Save Results in a CSV and .mat file
if sav
    t = tic;
    fprintf('Saving output in .mat and .csv files...');
    
    %
    [FDIRS , FNAMES , ~] = cellfun(@(x) fileparts(x), FPATHS, 'UniformOutput', 0);
    
    % Extract filename information
    ids  = {'UID' , 'Genotype'};
    vals = cellfun(@(x) getNameID(FNAMES, x), ids, 'UniformOutput', 0);
    
    % Extract PC scores
    allPCA = [ps , pt , ps , pr];
    scrs   = arrayfun(@(x) x.PCAScores, allPCA, 'UniformOutput', 0);
    
    % Store data in structure and table 
    flds = {'UID' , 'Genotype' , 'ShoulderPCs' , 'TipPCs' , ...
        'WidthPCs' , 'NormalizedPCs'};
    strc = cell2struct([vals , scrs] , flds, 2);
    tbl  = struct2table(strc);
    
    %
    dout = sprintf('%s_CarrotPCA_%dCarrots_%dGenotypes', ...
        tdate, numel(DSTS), numel(FDIRS));
    tnm1 = sprintf('%s/%s.csv', rootDir, dout);
    writetable(tbl, tnm1, 'FileType', 'text');
    
    tnm2 = sprintf('%s/%s', rootDir, dout);
    writetable(tbl, tnm2, 'FileType', 'spreadsheet');
    
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


