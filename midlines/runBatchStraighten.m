%%

% parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/model-validation/2019-groundtruth/original';
parentDir = '/Volumes/backup/images/gwas/pi-2019/rep_1/processed/raw';

% Set variables for extractor
vis = 0;
saveData = 0;
saveFigs = 0;

% Set variables for wrapper
straightenedMasks = 0; % smsk
widthProfile = 1; % dsts = 'vector' or barplot = 'barplot'
widthTicks = 0; % normals, not 
midlineOverlay = 1; % pmsk with mline and tcrd

% How far into the directory tree do you want to go
depth = 1

%%

batchStraighten(parentDir, vis, saveData, saveFigs, straightenedMasks, widthProfile, widthTicks, midlineOverlay, depth)

%%
