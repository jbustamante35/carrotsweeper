function [T , K , P , Q , error_files] = curvatureAnalysis(FPATHS, save_data, save_imgs, par, splt, num_pcs, out_pct, out_dim, fidx, vis)
%% curvatureAnalysis: run full pipeline to compute and plot curvatures
%
% Usage:
%   [T , K , P , Q , error_files] = ...
%       curvatureAnalysis(FPATHS, save_data, save_imgs, par, splt, ...
%       num_pcs, out_pct, out_dim, fidx, vis)
%
% Input:
%   FPATHS: cell array of file paths to images
%   save_data: boolean to save curvature values and PCA into csv file
%   save_imgs: boolean to save figures as png images
%   par: run silently [0], with parallel processing [1], or with verbosity [2]
%   splt: split by upper and lower sections (default: 0)
%   num_pcs:
%   out_pct: percentage of outliers to omit from PCA
%   out_dim:
%   fidx: figure handle index to plot data onto (default: 1)
%   vis: visualize results or skip (can't visualize if running with parallel)
%
% Output:
%   T: table containing curvature values as strings of each region and section
%   K: structure containing curvature values in numerical format
%   P: full results from PCA
%   Q: PCA results with outliers omitted
%   error_files: cell array of file names that caused errors
%

%% Set default parameters and options
if nargin < 2;  save_data = 0;        end
if nargin < 3;  save_imgs = 0;        end
if nargin < 4;  par       = 0;        end
if nargin < 5;  splt      = 0;        end
if nargin < 6;  num_pcs   = 5;        end
if nargin < 7;  out_pct   = [5 , 95]; end
if nargin < 8;  out_dim   = 1;        end
if nargin < 9;  fidx      = 1;        end
if nargin < 10; vis       = 0;        end

% Default parameters
EXCLUDE_COLUMNS = 15; % Number of left-most columns to exclude
SMOOTH_FACTOR   = 10; % Curvature smoothing factor
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
        ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));

        % NOTE: Can't visualize result in parallel mode!
        parfor f = 1 : numPaths
            try
                tAll = tic;
                fp   = FPATHS{f};
                nm   = sprintf('%s | %s', GENOTYPE{f}, snms{f});
                fprintf('%s\n%s:\n%s\n', sprA, nm, sprB);

                % Extract image from file path
                t = tic;
                str  = sprintf('Extracting Image');
                msks = imcomplement(double(logical(imread(fp))));
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));

                % Compute curvatures
                t   = tic;
                str = sprintf('Computing Curvatures...');
                sk{f} = computeCurvatures(msks, splt, EXCLUDE_COLUMNS, ...
                    SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP);
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));

                str = sprintf('Finished %04d of %04d', f, numPaths);
                fprintf('%s\n%s%s [ %.02f sec ]\n%s] \n', ...
                    sprB, str, ellp(str), toc(tAll), sprA);

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
        ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));

        for f = 1 : numPaths
            try
                tAll = tic;
                fp   = FPATHS{f};
                nm   = sprintf('%s | %s', GENOTYPE{f}, snms{f});
                fprintf('%s\n%s:\n%s\n', sprA, nm, sprB);

                % Extract image from file path
                t = tic;
                str  = sprintf('Extracting Image');
                msks = imcomplement(double(logical(imread(fp))));
                fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));

                % Compute curvatures
                t   = tic;
                str = sprintf('Computing Curvatures...');
                [sk{f} , sc, sm, sskel] = ...
                    computeCurvatures(msks, splt, EXCLUDE_COLUMNS, ...
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
                    sprB, str, ellp(str), toc(tAll), sprA);
            catch err
                % Store filename of error and continue loop
                error_files{error_count} = ...
                    sprintf('%s_%s', GENOTYPE{f}, snms{f});
                fprintf(2, '\n\n%s\n\nError with %s | Total Errors: %d\n', ...
                    err.getReport, error_files{error_count}, error_count);

                error_count = error_count + 1;
                continue;
            end
        end

    otherwise
        %% Run silently using cell functions
        msks                    = cellfun(@(x) imcomplement(double(logical(imread(x)))), ...
            FPATHS, 'UniformOutput', 0);
        [sk , sc , sm , sskel]  = cellfun(@(x) computeCurvatures(x, ...
            splt, EXCLUDE_COLUMNS, SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP), ...
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
K                   = cat(1, sk{:});
[P , Q , s , t , w] = curvaturePCA(K, splt, num_pcs, out_pct, out_dim);

%% Store curvatures and curvature figures into individual folders
% Iterate through regions and sections for straightened and binary curvatures
sK = cellfun(@(k,n) processCurvatures(k,n, splt), sk, snms, 'UniformOutput', 0);

% Convert to tables
S = cat(1, sK{:});
T = struct2table(S);

%% Add PCA Scores and Miscllaneous Curvature Data to table T
% Split upper-lower sections by left-right, then add to table
if isempty(Q)
    scrss = P.shoulder.PCAScores;
    scrst = P.tip.PCAScores;
    scrsw = P.whole.PCAScores;

    % Eigenvectors and Means
    evecs = P.shoulder.EigVecs;
    evect = P.tip.EigVecs;
    evecw = P.whole.EigVecs;
    mnss  = P.shoulder.MeanVals;
    mnst  = P.tip.MeanVals;
    mnsw  = P.whole.MeanVals;

else
    scrss = Q.shoulder.PCAScores;
    scrst = Q.tip.PCAScores;
    scrsw = Q.whole.PCAScores;

    % Eigenvectors and Means
    evecs = Q.shoulder.EigVecs;
    evect = Q.tip.EigVecs;
    evecw = Q.whole.EigVecs;
    mnss  = Q.shoulder.MeanVals;
    mnst  = Q.tip.MeanVals;
    mnsw  = Q.whole.MeanVals;
end

if splt
    % Sections are split upper and lower and need to be reshaped from a
    % 2-dimensional matrix to a 3-dimensional matrix
    %
    % [ u11 u12 ... u1p ]    [ u11 l11 ]   [ u12 l12 ]  ...  [ u1p l1p ]
    % [ l11 l12 ... l1p ]    [ u21 l21 ]   [ u22 l22 ]  ...  [ u2p l2p ]
    % [ u21 u22 ... u2p ]        ...          ...       ...      ...
    % [ l21 l22 ... l2p ] => [ un1 ln1 ]   [ un2 ln2 ]  ...  [ unp lnp ]
    %   ...     ...             PC 1          PC 2      ...     PC p
    % [ un1 un2 ... unp ]
    % [ ln1 ln2 ... lnp ]

    % Shoulders
    [cnvps , sumks , avgks] = reshapeSplitCurvatures(scrss, s);

    T.shoulder_upper_sum = sumks(:,1);
    T.shoulder_upper_avg = avgks(:,1);
    T.shoulder_lower_sum = sumks(:,2);
    T.shoulder_lower_avg = avgks(:,2);

    % Tips
    [cnvpt , sumkt , avgkt] = reshapeSplitCurvatures(scrst, t);

    T.tip_upper_sum = sumkt(:,1);
    T.tip_upper_avg = avgkt(:,1);
    T.tip_lower_sum = sumkt(:,2);
    T.tip_lower_avg = avgkt(:,2);

    % Whole
    [cnvpw , sumkw , avgkw] = reshapeSplitCurvatures(scrsw, w);

    T.whole_upper_sum = sumkw(:,1);
    T.whole_upper_avg = avgkw(:,1);
    T.whole_lower_sum = sumkw(:,2);
    T.whole_lower_avg = avgkw(:,2);

    for i = 1 : num_pcs
        T.(sprintf('shoulder_upper_pc%d', i)) = squeeze(cnvps(:,1,i));
        T.(sprintf('shoulder_lower_pc%d', i)) = squeeze(cnvps(:,2,i));
        T.(sprintf('tip_upper_pc%d', i))      = squeeze(cnvpt(:,1,i));
        T.(sprintf('tip_lower_pc%d', i))      = squeeze(cnvpt(:,1,i));
        T.(sprintf('whole_upper_pc%d', i))    = squeeze(cnvpw(:,1,i));
        T.(sprintf('whole_lower_pc%d', i))    = squeeze(cnvpw(:,1,i));
    end

else
    % Shoulders
    T.shoulder_sum = sum(s,2);
    T.shoulder_avg = mean(s,2);

    % Tips
    T.tip_sum = sum(t,2);
    T.tip_avg = mean(t,2);

    % Whole
    sumkw = sum(w,2);
    avgkw = mean(w,2);

    T.whole_sum = sumkw(:,1);
    T.whole_avg = avgkw(:,1);

    for i = 1 : num_pcs
        T.(sprintf('shoulder_pc%d', i)) = scrss(:,i);
        T.(sprintf('tip_pc%d', i))      = scrst(:,i);
        T.(sprintf('whole_pc%d', i))    = scrsw(:,i);
    end
end

%% Ouput table into csv/xls file and files that caused errors in csv file
if save_data
    if ~isfolder(OUT_DIR); mkdir(OUT_DIR); end

    % Save Curvatures and PC Scores
    tnm = sprintf('%s%s%s_Curvatures_%04dCarrots', ...
        OUT_DIR, filesep, tdate, numel(S));
    writetable(T, [tnm , '.csv'], 'FileType', 'text');

    % Save eigenvectors and means
    estrs = sprintf('%s_Curvatures_ShoulderVectors', tdate);
    estrt = sprintf('%s_Curvatures_TipVectors', tdate);
    estrw = sprintf('%s_Curvatures_WholeVectors', tdate);

    enms = struct('EigVecs', evecs, 'Means', mnss');
    enmt = struct('EigVecs', evect, 'Means', mnst');
    enmw = struct('EigVecs', evecw, 'Means', mnsw');

    etbls = struct2table(enms);
    etblt = struct2table(enmt);
    etblw = struct2table(enmw);

    etnms = sprintf('%s/%s.csv', OUT_DIR, estrs);
    etnmt = sprintf('%s/%s.csv', OUT_DIR, estrt);
    etnmw = sprintf('%s/%s.csv', OUT_DIR, estrw);

    writetable(etbls, etnms, 'FileType', 'text');
    writetable(etblt, etnmt, 'FileType', 'text');
    writetable(etblw, etnmw, 'FileType', 'text');

    % Error files
    enm = sprintf('%s%s%s_errorfiles_%04dFiles', ...
        OUT_DIR, filesep, tdate, numel(error_files));
    E   = cell2table(error_files);
    writetable(E, [enm , '.csv'], 'FileType', 'text');

end
end

function [P , Q , s , t , w] = curvaturePCA(K, splt, num_pcs, out_pct, out_dim)
%% curvaturePCA: run PCA on curvatures and omit outliers
if splt
    % Upper and lower sections are spilt such that the upper section data are
    % the first half of the matrix, and the lower section is the second half
    s = [arrayfun(@(k) k.shoulder.upper', K, 'UniformOutput', 0) ; ...
        arrayfun(@(k) flipud(k.shoulder.lower)', K, 'UniformOutput', 0)];
    t = [arrayfun(@(k) k.tip.upper', K, 'UniformOutput', 0) ; ...
        arrayfun(@(k) flipud([k.tip.lower ; 0])', K, 'UniformOutput', 0)];
    w = [arrayfun(@(k) k.whole.upper', K, 'UniformOutput', 0) ; ...
        arrayfun(@(k) flipud(k.whole.lower)', K, 'UniformOutput', 0)];

else
    % Upper and Lower sections are unsplit and stitched together
    s = arrayfun(@(k) k.shoulder', K, 'UniformOutput', 0);
    t = arrayfun(@(k) k.tip', K, 'UniformOutput', 0);
    w = arrayfun(@(k) k.whole', K, 'UniformOutput', 0);
end

s = cat(1, s{:});
t = cat(1, t{:});
w = cat(1, w{:});

% Run PCA and omit outliers [store scores before omitting outliers]
P.shoulder = pcaAnalysis(s, num_pcs, 0, 'shoulders_split');
P.tip      = pcaAnalysis(t, num_pcs, 0, 'tip_split');
P.whole    = pcaAnalysis(w, num_pcs, 0, 'whole_split');

% Remove Top and Bottom Outliers
if ~isempty(out_pct)
    Q.shoulder = pcaOmitOutliers(P.shoulder, out_pct, out_dim);
    Q.tip      = pcaOmitOutliers(P.tip,      out_pct, out_dim);
    Q.whole    = pcaOmitOutliers(P.whole,    out_pct, out_dim);
else
    Q = [];
end

end

function [cnvp , sumk , avgk] = reshapeSplitCurvatures(scrs, k)
%% reshapeCurvatures: reshape split curvature data from 2D to 3D
cnvp = reshape(scrs,     [size(scrs,1) / 2 , 2 , size(scrs,2)]);
sumk = reshape(sum(k,2),  [size(k,1) / 2 , 2]);
avgk = reshape(mean(k,2), [size(k,1) / 2 , 2]);

end
