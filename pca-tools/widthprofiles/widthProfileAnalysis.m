function [pw, ps, pt, pr] =  widthProfileAnalysis(rootDir, maskDir, sav)
%% widthProfileAnalysis: mutliple PCA analyses of width profiles
% Get PCA on Widths, Length-Width Normalized Profile, Tips, Shoulders.
% A neat pipeline to do all this shit from raw straightened masks. Output is in
% the form of PCA objects to all of these, and other data too but I'm not sure
% what else to store that can't be extrapolated from the PCA data.
%
% Usage:
%   [tips, shldr, wids] =  widthProfileAnalysis(rootDir, maskDir, sav)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Raw Width Profiles from Straightened Masks
[DSTS , ~, FPATHS] = loadWidthProfiles(rootDir, maskDir);

%% Prepare Normalized Width Profiles and Split Tips and Shoulders
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

%% PCA on Width, Tips, and Shoulders
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

%% Graham-Schmidt Analyses to perform orthonormalization of width profiles
% Covariance in data using Canonical Correlations Analysis
COVARMETHOD = 'cca';
pscrs       = pw.PCAScores;
wids_lens   = [wids , lens];
[bb, preB]  = pcaRegression(wids_lens, pscrs, COVARMETHOD);

% Use regressors to remove length and width information from dataset
H = grahamSchmidt(bb(:,1), bb(:,2), pscrs);

%% PCA on length-width normalized profiles
wrnm = sprintf('%dNormalizedLengthWidth', numel(DSTS));
pr   = pcaAnalysis(H, pcrw, sav, wrnm, 0);

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

end

function vals = getNameID(FNAMES, id)
%% getNameID: extract values from an id in a filename
%
expr = sprintf('%s_(?<id>.*?)}', id);
val  = regexpi(FNAMES, expr, 'names');
vals = cellfun(@(x) char(x.id), val, 'UniformOutput', 0);

end


