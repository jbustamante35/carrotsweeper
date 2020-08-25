function [T , K , P , Q , error_files] = curvatureAnalysis(FPATHS, save_data, save_imgs, par, num_pcs, out_pct, out_dim, fidx, vis)
%% curvatureAnalysis: run full pipeline to compute and plot curvatures
%
% Usage:
%   [T , K , P , Q , error_files] = curvatureAnalysis(FPATHS, save_data, save_imgs, ...
%            par, num_pcs, out_pct, out_dim, fidx, vis)
%
% Input:
%   FPATHS: cell array of file paths to images
%   save_data: boolean to save curvature values and PCA into csv file
%   save_imgs: boolean to save figures as png images
%   par: run silently [0], with parallel processing [1], or with verbosity [2]
%   out_pct: percentage of outliers to omit from PCA
%   fidx: figure handle index to plot data onto (default: 1)
%   vis: visualize results or skip (can't visualize if running with parallel)
%
% Output:
%   T: table containing curvature values of each region and section
%   K: structure containing curvature values
%   P: full results from PCA
%   Q: PCA results with outliers omitted
%   error_files: cell array of file names that caused errors
%

%% Set default parameters and options
if nargin < 2
    save_data = 0;
    save_imgs = 0;
    fidx      = 1;
    vis       = 0;
    par       = 0;
    num_pcs   = 5;
    out_pct   = [5 , 95];
    out_dim   = 1;
end

% Default parameters
SPLIT_REGIONS   = 1;  % Split by upper and lower sections [1]
EXCLUDE_COLUMNS = 15; % Number of left-most columns to exclude
SMOOTH_FACTOR   = 14; % Curvature smoothing factor
SIZE_SHOULDER   = 50; % Number of coordinates to sample for shoulders
SIZE_TIP        = 50; % Number of coordinates to sample for tips
DSTR            = 0;  % Show distribution of curvatures
K_DIR_NAME      = 'curvature-mask';
OUT_DIR_NAME    = 'Output';

% Details for figure
[kdirs , snms] = cellfun(@(x) fileparts(x), FPATHS, 'UniformOutput', 0);
K_DIR          = cellfun(@(d) [fileparts(d) , filesep , K_DIR_NAME], ...
    kdirs, 'UniformOutput', 0);
[~ , GENOTYPE] = cellfun(@(d) fileparts(fileparts(d)), ...
    kdirs, 'UniformOutput', 0);

% Output directory
OUT_DIR = [fileparts(fileparts((K_DIR{1}))) , filesep , OUT_DIR_NAME];

%% Curvature analysis
numPaths    = numel(FPATHS);
error_files = cell(1, numPaths);
error_count = 1;

% Get masks from file paths and compute curvatures
switch par
    case 1
        %% Run with parallel processing (be careful with RAM!)
        fprintf('NOTE: Can''t visualize results in parallel mode!');
        
        halfCores = ceil(feature('numcores') / 2);
        setupParpool(halfCores, 0);
        
        sk   = cell(numPaths, 1);
        sepA = repmat('=', [1 , 80]);
        sepB = repmat('-', [1 , 80]);
        ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));
        
        % NOTE: Can't visualize result in parallel mode!
        parfor f = 1 : numPaths
            try
                tAll = tic;
                fp   = FPATHS{f};
                nm   = sprintf('%s | %s', GENOTYPE{f}, snms{f});
                fprintf('%s\n%s:\n%s\n', sepA, nm, sepB);
                
                % Extract image from file path
                t = tic;
                str  = sprintf('Extracting Image');
                msks = imcomplement(double(logical(imread(fp))));
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
                
                % Compute curvatures
                t   = tic;
                str = sprintf('Computing Curvatures...');
                sk{f} = computeCurvatures(msks, SPLIT_REGIONS, EXCLUDE_COLUMNS, ...
                    SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP);
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
                
                str = sprintf('Finished %04d of %04d', f, numPaths);
                fprintf('%s\n%s%s [ %.02f sec ]\n%s] \n', ...
                    sepB, str, ellp(str), toc(tAll), sepA);
                
            catch
                % Store filename of error and continue loop
                error_files{f} = sprintf('%s_%s', GENOTYPE{f}, snms{f});
                fprintf(2, 'Error with %s\n', error_files{f});
                continue;
            end
        end
    case 2
        %% Run in for loop with verbose output
        sk   = cell(numPaths, 1);
        sepA = repmat('=', [1 , 80]);
        sepB = repmat('-', [1 , 80]);
        ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));
        
        for f = 1 : numPaths
            try
                tAll = tic;
                fp   = FPATHS{f};
                nm   = sprintf('%s | %s', GENOTYPE{f}, snms{f});
                fprintf('%s\n%s:\n%s\n', sepA, nm, sepB);
                
                % Extract image from file path
                t = tic;
                str  = sprintf('Extracting Image');
                msks = imcomplement(double(logical(imread(fp))));
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
                
                % Compute curvatures
                t   = tic;
                str = sprintf('Computing Curvatures...');
                [sk{f} , sc, sm, sskel] = ...
                    computeCurvatures(msks, SPLIT_REGIONS, EXCLUDE_COLUMNS, ...
                    SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP);
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
                
                % Plot curvatures and save image
                if vis
                    t   = tic;
                    str = sprintf('Plotting Curvatures...');
                    plotCurvature(sskel, sk{f}, sc, sm, ...
                        SMOOTH_FACTOR, EXCLUDE_COLUMNS, save_imgs, snms{f}, K_DIR{f}, fidx, DSTR);
                    fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
                    
                end
                
                str = sprintf('Finished %04d of %04d', f, numPaths);
                fprintf('%s\n%s%s [ %.02f sec ]\n%s] \n', ...
                    sepB, str, ellp(str), toc(tAll), sepA);
            catch
                % Store filename of error and continue loop
                error_files{error_count} = sprintf('%s_%s', GENOTYPE{f}, snms{f});
                fprintf(2, 'Error with %s | Total Errors: %d\n', ...
                    error_files{error_count}, error_count);
                
                error_count = error_count + 1;
                continue;
            end
        end
        
    otherwise
        %% Run silently using cell functions
        msks                    = cellfun(@(x) imcomplement(double(logical(imread(x)))), ...
            FPATHS, 'UniformOutput', 0);
        [sk , sc , sm , sskel]  = cellfun(@(x) computeCurvatures(x, ...
            SPLIT_REGIONS, EXCLUDE_COLUMNS, SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP), ...
            msks, 'UniformOutput', 0);
        
        if vis
            % Plot result
            cellfun(@(s,k,c,m,nm,do) plotCurvature(s, k, c, m, ...
                SMOOTH_FACTOR, EXCLUDE_COLUMNS, save_imgs, nm, do, fidx), ...
                sskel, sk, sc, sm, snms, K_DIR, 'UniformOutput', 0);
        end
end

%% Concatenate error files and remove from dataset
error_files = error_files(~cellfun(@isempty, error_files));
echk        = cellfun(@(e) find(cell2mat(cellfun(@(n) contains(e, n), ...
    snms, 'UniformOutput', 0))), error_files, 'UniformOutput', 0);
echk        = cat(1, echk{:});
snms(echk)  = [];
sk(echk)    = [];

%% Run PCA on shoulder and tip curvatures
% Combine upper and lower into same dataset
K       = cat(1, sk{:});
[P , Q] = curvaturePCA(K, num_pcs, out_pct, out_dim);

%% Store curvatures and curvature figures into individual folders
% Iterate through regions and sections for straightened and binary curvatures
sK = cellfun(@(k,n) processCurvatures(k,n), sk, snms, 'UniformOutput', 0);

% Convert to tables
S = cat(1, sK{:});
T = struct2table(S);

%% Add PCA Scores and Miscllaneous Curvature Data to table T
% Split upper-lower sections by left-right, then add to table
scrss = P.shoulder.PCAScores;
scrst = P.tip.PCAScores;

% Shoulders
cnvps = im2colF(scrss(:,1), [2 , 1], [2 , 1])';
sumks = im2colF(sum(s,2), [2 , 1], [2 , 1])';
avgks = im2colF(mean(s,2), [2 , 1], [2 , 1])';

T.shoulder_upper_pcs = cnvps(:,1);
T.shoulder_upper_sum = sumks(:,1);
T.shoulder_upper_avg = avgks(:,1);
T.shoulder_lower_pcs = cnvps(:,2);
T.shoulder_lower_sum = sumks(:,2);
T.shoulder_lower_avg = avgks(:,2);

% Tips
cnvpt = im2colF(scrst(:,1), [2 , 1], [2 , 1])';
sumkt = im2colF(sum(t,2), [2 , 1], [2 , 1])';
avgkt = im2colF(mean(t,2), [2 , 1], [2 , 1])';

T.tip_upper_pcs = cnvpt(:,1);
T.tip_upper_sum = sumkt(:,1);
T.tip_upper_avg = avgkt(:,1);
T.tip_lower_pcs = cnvpt(:,2);
T.tip_lower_sum = sumkt(:,2);
T.tip_lower_avg = avgkt(:,2);

%% Ouput table into xls file
if save_data
    if ~isfolder(OUT_DIR)
        mkdir(OUT_DIR);
    end
    
    tnm = sprintf('%s%s%s_Curvatures_%04dCarrots', ...
        OUT_DIR, filesep, tdate, numel(S));
    writetable(T, [tnm , '.csv'], 'FileType', 'text');
    writetable(T, tnm, 'FileType', 'spreadsheet');
end
end

function [P , Q] = curvaturePCA(K, num_pcs, out_pct, out_dim)
%% curvaturePCA: run PCA on curvatures and omit outliers
s = arrayfun(@(k) [k.shoulder.upper , flipud(k.shoulder.lower)]', ...
    K, 'UniformOutput', 0);
t = arrayfun(@(k) [k.tip.upper , flipud([k.tip.lower ; 0])]', ...
    K, 'UniformOutput', 0);

s = cat(1, s{:});
t = cat(1, t{:});

% Run PCA and omit outliers [store scores before omitting outliers]
P.shoulder = pcaAnalysis(s, num_pcs, 0, 'shoulders_split', 0);
P.tip      = pcaAnalysis(t, num_pcs, 0, 'tip_split', 0);

% Remove Top and Bottom Outliers
Q.shoulder = pcaOmitOutliers(P.shoulder, out_pct, out_dim);
Q.tip      = pcaOmitOutliers(P.tip, out_pct, out_dim);

end
