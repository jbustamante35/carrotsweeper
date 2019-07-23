%%

% parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/model-validation/2019-groundtruth/original';
% parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated';
parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6166';

% Set variables for extractor
vis = 0;
saveData = 0;
saveFigs = 0;

% Set variables for wrapper
straightenedMasks = 1; % smsk
widthProfile = 0; % dsts = 'vector' or barplot = 'barplot'
widthTicks = 0; % normals, not 
midlineOverlay = 1; % pmsk with mline and tcrd

% How far into the directory tree do you want to go
depth = 76

%%

batchStraighten(parentDir, vis, saveData, saveFigs, straightenedMasks, widthProfile, widthTicks, midlineOverlay, depth)

%%


%% 
binaryMaskPath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6619/binary-masks/';
genoPath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6619/';
