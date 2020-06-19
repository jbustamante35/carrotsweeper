function W = prepWidthAnalysis(varargin)
%% prepWidthAnalysis:
% Description
%
% Usage:
%    W = prepWidthAnalysis(varargin)
%
% Input:
%    rootDir: path to root directory of image sub-directories
%    maskDir: name of the sub-directory containing straightened-masks
%    varargin: instructions on how to prepare width profiles
%
% Output:
%   W: structure containing width profiles and regressors, and save filename
%
% Properties in output structure W:
%   Name: filename to save as
%   Profiles: width profiles prepared for PCA
%   WidthRegressor: width regressors to convert from normalized to raw value
%   LengthRegressor: length regressors  to convert from normalized to raw value
%
% Author Julian Bustamante <jbustamante@wisc.edu>
%

%% Parse inputs to obtain property values
prps = {'RootDir' , 0 ; 'MaskDir' , 'straight-masks' ; 'Profiles', [] ; ...
    'LengthNormalize' , 0 ; 'WidthNormalize' , 0 ; ...
    'LengthInterp' , 1000 ; 'GrahamSchmidt' , 0 ; ...
    'ClipSection' , 'whole' ; 'ClipSize' , 150 ; 'Save', 0};
deflts = prps(~cellfun(@isempty, prps(:,2)),:);
vargs  = {varargin};
args   = parseConstructorInput(prps(:,1), deflts, vargs);

%% Get Raw Width Profiles from Straightened Masks or load Profiles
if args.RootDir
    t = tic;
    fprintf('Loading width profiles from %s/%s...', rootDir, maskDir);
    
    [DSTS , ~, ~] = loadWidthProfiles(rootDir, maskDir);
    ttlWids       = numel(DSTS);
    
    fprintf('DONE! [%.02f sec]\n', toc(t));
    
elseif ~isempty(args.Profiles)
    DSTS    = args.Profiles;
    ttlWids = numel(DSTS);
end

%% Prepare Normalized/Non-Normalized Width Profiles of Full, Tips, or Shoulders
t = tic;
fprintf('Prepping %d profiles [Section = %s|LengthNorm = %d|WidthNorm = %d|GM = %d]...', ...
    ttlWids, args.ClipSection, args.LengthNormalize, args.WidthNormalize, ...
    args.GrahamSchmidt);

%% Determine length and width normalization methods
nrmL       = args.LengthNormalize;
nrmW       = args.WidthNormalize;
clipinterp = args.LengthInterp;

if nrmL
    ntL = 'Normalized';
else
    ntL = 'Original';
end

if nrmW
    ntW = 'Normalized';
else
    ntW = 'Original';
end

%% Use full width profile or clip off shoulders or tips
sectiontype   = args.ClipSection;
clip          = 1 : args.ClipSize;
[wids , lens] = deal(0);

switch sectiontype
    case 'whole'
        %% Analyze full width profile
        [P , wids , lens] = prepWidthProfiles(DSTS, nrmL, nrmW, clipinterp);
        
        % Graham-Schmidt Analyses to perform orthonormalization of width profiles
        if args.GrahamSchmidt
            tt = tic;
            fprintf('Orthonormalization on %d width profiles...', ttlWids);
            
            % Graham-Shmidt [width and length normalization]
            sectiontype = 'GrahamSchmidt';
            COVARMETHOD = 'cca';
            gmpcs       = args.GrahamSchmidt;
            
            pw          = myPCA(P, gmpcs);
            scrs        = pw.PCAScores;
            wids_lens   = [wids , lens];
            [bb, ~]     = pcaRegression(wids_lens, scrs, COVARMETHOD);
            
            % Use regressors to remove length and width information from dataset
            P = grahamSchmidt(bb(:,1), bb(:,2), scrs);
            
            fprintf('DONE! [%.02f sec]...', toc(tt));
        end
        
    case 'shoulders'
        %% Clip off shoulders (first x section of width profile)
        FLPS = cellfun(@fliplr, DSTS, 'UniformOutput', 0);
        SHLD = cellfun(@(x) x(clip), FLPS, 'UniformOutput', 0);
        P    = fliplr(prepWidthProfiles(SHLD, nrmL, nrmW, clipinterp));
        
    case 'tips'
        %% Clip off tips (last x section of width profile)
        TIPS = cellfun(@(x) x(clip), DSTS, 'UniformOutput', 0);
        P    = prepWidthProfiles(TIPS, nrmL, nrmW, clipinterp);
        
    otherwise
        fprintf(2, 'Error with ClipSection %s [full|tips|shoulders]\n', ...
            args.ClipSection);
        return;
end

fprintf('DONE! [%.02f sec]\n', toc(t));

%% Output structure
nm = sprintf('%s_%dProfiles_%sLengths_%sWidths_%s', ...
    tdate, ttlWids, ntL, ntW, sectiontype);
W  = struct('Name', nm, 'Profiles', P, 'LengthNorm', ntL, 'WidthNorm', ntW, ...
    'WidthRegressor', wids, 'LengthRegressor', lens);

if args.Save
    save(nm, '-v7.3', 'W');
end

end

function args = parseConstructorInput(prps, deflts, vargs)
%% Parse input parameters for Constructor method
p = inputParser;

% Replace empty argument with default value
emptyprps = cell(numel(prps), 1);
if ~isempty(deflts)
    matchIdxs            = cell2mat(cellfun(@(x) find(strcmp(prps, x)), ...
        deflts(:,1), 'UniformOutput', 0));
    emptyprps(matchIdxs) = deflts(:,2);
end

% Add all properties as empty
cellfun(@(x,d) p.addOptional(x, d), prps, emptyprps, 'UniformOutput', 0);

% Parse arguments and output into structure
p.parse(vargs{1}{:});
args = p.Results;
end
