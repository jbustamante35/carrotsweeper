function T = curvatureAnalysis(FPATHS, save_data, save_fig, fidx, par)
%% curvatureAnalysis: run full pipeline to compute and plot curvatures
%
% Usage:
%   T = curvatureAnalysis(SPATHS, save_data, save_fig, fidx, par)
%
% Input:
%   FPATHS: cell array of file paths to images
%   save_data: boolean to save curvatures into xls file
%   save_fig: boolean to save figures as png images
%   fidx: figure handle index to plot data onto (default: 1)
%   par: run silently [0], with parallel processing [1], or with verbosity [2]
%
% Output:
%   T: table containing curvature values of each region and section
%

%% Set default parameters and options
if nargin < 2
    save_data = 0;
    save_fig  = 0;
    fidx      = 1;
    par       = 0;
end

% Default parameters
SPLIT_REGIONS   = 1;
EXCLUDE_COLUMNS = 15;
SMOOTH_FACTOR   = 14;
SIZE_SHOULDER   = 50;
SIZE_TIP        = 50;
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
% Get masks from file paths and compute curvatures
switch par
    case 1
        %% Run with parallel processing (be careful with RAM!)
        fprintf(2, '\nDon''t use with parallelization yet\n\n');
        T = [];
        return;
        %     halfCores = ceil(feature('numcores') / 2);
        %     setupParpool(halfCores, 0);
        %
        %         numPaths = numel(FPATHS);
        %         sk       = cell(numPaths, 1);
        %
        %         sepA = repmat('=', [1 , 80]);
        %         sepB = repmat('-', [1 , 80]);
        %         ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));
        %         parfor f = 1 : numPaths
        %             tAll = tic;
        %             fp   = FPATHS{f};
        %             nm   = sprintf('%s | %s', GENOTYPE{f}, snms{f});
        %             fprintf('%s\n%s:\n%s\n', sepA, nm, sepB);
        %
        %             % Extract image from file path
        %             t = tic;
        %             str  = sprintf('Extracting Image');
        %             msks = imcomplement(double(logical(imread(fp))));
        %             fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
        %
        %             % Compute curvatures
        %             t   = tic;
        %             str = sprintf('Computing Curvatures...');
        %             [sk{f} , sc, sm, sskel] = ...
        %                 computeCurvatures(msks, SPLIT_REGIONS, EXCLUDE_COLUMNS, ...
        %                 SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP);
        %             fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
        %
        %             % Plot curvatures and save image
        %             t   = tic;
        %             str = sprintf('Plotting Curvatures...');
        %             plotCurvature(sskel, sk{f}, sc, sm, ...
        %                 SMOOTH_FACTOR, EXCLUDE_COLUMNS, save_fig, snms{f}, K_DIR{f}, fidx);
        %             fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
        %
        %             str = sprintf('Finished %04d of %04d', f, numPaths);
        %             fprintf('%s\n%s%s [ %.02f sec ]\n%s] \n', ...
        %                 sepB, str, ellp(str), toc(tAll), sepA);
        %         end
        
    case 2
        %% Run in for loop with verbose output
        numPaths = numel(FPATHS);
        sk       = cell(numPaths, 1);
        
        sepA = repmat('=', [1 , 80]);
        sepB = repmat('-', [1 , 80]);
        ellp = @(x) repmat('.', 1, 80 - (length(x) + 13));
        for f = 1 : numPaths
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
            t   = tic;
            str = sprintf('Plotting Curvatures...');
            plotCurvature(sskel, sk{f}, sc, sm, SMOOTH_FACTOR, EXCLUDE_COLUMNS, ...
                save_fig, snms{f}, K_DIR{f}, fidx, 0);
            fprintf('%s%s [ %.02f sec ] \n', str, ellp(str), toc(t));
            
            str = sprintf('Finished %04d of %04d', f, numPaths);
            fprintf('%s\n%s%s [ %.02f sec ]\n%s] \n', ...
                sepB, str, ellp(str), toc(tAll), sepA);
        end
        
    otherwise
        %% Run silently using cell functions
        msks                    = cellfun(@(x) imcomplement(double(logical(imread(x)))), ...
            FPATHS, 'UniformOutput', 0);
        [sk , sc , sm , sskel]  = cellfun(@(x) computeCurvatures(x, ...
            SPLIT_REGIONS, EXCLUDE_COLUMNS, SMOOTH_FACTOR, SIZE_SHOULDER, SIZE_TIP), ...
            msks, 'UniformOutput', 0);
        
        % Plot result
        cellfun(@(s,k,c,m,nm,do) plotCurvature(s, k, c, m, SMOOTH_FACTOR, ...
            EXCLUDE_COLUMNS, save_fig, nm, do, fidx, 0), ...
            sskel, sk, sc, sm, snms, K_DIR, 'UniformOutput', 0);
end

%% Store curvatures and curvature figures into individual folders
% Iterate through regions and sections for straightened and binary curvatures
sK  = cellfun(@(k,n) processCurvatures(k,n), sk, snms, 'UniformOutput', 0);

% Convert to table
S = cat(1, sK{:});
T = struct2table(S);

%% Ouput table into xls file
if save_data
    if ~isfolder(OUT_DIR)
        mkdir(OUT_DIR);
    end
    
    tnm = sprintf('%s%s%s_Curvatures_%04dCarrots', ...
        OUT_DIR, filesep, tdate, numel(S));
    writetable(T, tnm, 'FileType', 'spreadsheet');
end
end
