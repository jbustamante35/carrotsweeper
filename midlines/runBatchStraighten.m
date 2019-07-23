%%

% parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/model-validation/2019-groundtruth/original';
<<<<<<< HEAD
% parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated';
parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6166';
=======
parentDir = '/Volumes/backup/images/gwas/pi-2019/rep_1/processed/raw';
>>>>>>> c32c7bff76a4cd10d54ade0fc5cf6523e045230c

% Set variables for extractor
vis = 0;
saveData = 0;
saveFigs = 0;

% Set variables for wrapper
<<<<<<< HEAD
straightenedMasks = 1; % smsk
widthProfile = 0; % dsts = 'vector' or barplot = 'barplot'
=======
straightenedMasks = 0; % smsk
widthProfile = 1; % dsts = 'vector' or barplot = 'barplot'
>>>>>>> c32c7bff76a4cd10d54ade0fc5cf6523e045230c
widthTicks = 0; % normals, not 
midlineOverlay = 1; % pmsk with mline and tcrd

% How far into the directory tree do you want to go
<<<<<<< HEAD
depth = 76
=======
depth = 1
>>>>>>> c32c7bff76a4cd10d54ade0fc5cf6523e045230c

%%

batchStraighten(parentDir, vis, saveData, saveFigs, straightenedMasks, widthProfile, widthTicks, midlineOverlay, depth)

%%
<<<<<<< HEAD


%% 
binaryMaskPath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6619/binary-masks/';
genoPath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/gwas-curated/NSL-6619/';
=======
>>>>>>> c32c7bff76a4cd10d54ade0fc5cf6523e045230c
