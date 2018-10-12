function midlines = carrotMidlinesPipeline(image_directory, reset_directory)
%% carrotMidlines: format midline data from phytoMorph algorithm
% I'm not sure what this will do just yet. Something along the lines of:
% _/ Crop binary mask image to push shoulder to edge of image
% _/ Copy processed mask 5 times to make a "time-lapse"
% X Run customized phytoMorph midline algorithm [running from CyVerse is overkill]
% X Extract midline data from phytoMorph output
%
% Usage:
%   midlines = carrotMidlinesPipeline(image_directory, reset_directory)
%
% Input:
%   image_directory: directory containing binary mask images to compute midlines
%   reset_directory: boolean to run through pipeline normally (0) or reset directory (1)
%
% Output:
%   midlines: outputted midline data from phytoMorph algorithm
%

%% Process images and Set-up directory structure
% Create image storage object and create backup directory
O = imageDatastore(image_directory);
workDir = sprintf('%s_midlines', image_directory);
mkdir(workDir);
copyfile(image_directory, workDir);
cd(workDir);
I = imageDatastore(workDir);

% Crop out edges of bw masks
img = arrayfun(@(x) push2edge(I.readimage(x)), 1 : numel(I.Files), 'UniformOutput', 0);
for m = 1 : numel(I.Files)
    replaceOrignalWithCrop(I.Files{m}, img{m});
end

% Create directories for each image and rename to 1.ext [ext == image extension]
[dirNames, ext] = cellfun(@(x) formatDirectory(x), I.Files, 'UniformOutput', 0);

% Go into each new directory and copy file n times
cd(workDir);
num_copies = 5;
cellfun(@(x) copyImage(workDir, ext{1}, x, num_copies), dirNames, 'UniformOutput', 0);

%% RESET OPTION IF SOMETHING SCREWS UP
if reset_directory
    cellfun(@(x) resetDir(workDir, x), dirNames, 'UniformOutput', 0);
    for m = 1 : numel(O.Files)
        replaceOrignalWithCrop(I.Files{m}, O.readimage(m));
    end
end

%% Run images through phytoMorph algorithm and process output data [TODO]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Customize phytoMorph algorithm to work locally
%                           Extract midline data from phytoMorph output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M        = cellfun(@(x) phytoMorphMidline(workDir, x), dirNames, 'UniformOutput', 0);
midlines = cellfun(@(x) extractMidlineData(x), M, 'UniformOutput', 0);
cd('..');
end

function [dirName, ext] = formatDirectory(fin)
%% formatDirectory: sub-function to create new directory from filename
% Extract filename of image defined by fin
% Create directory named fin and move image into new directory
% Change name of image to 1.ext [ext == image extension]
fName   = getDirName(fin);
ext     = fName((strfind(fName, '.') + 1) : end);
dirName = fName(1 : (strfind(fName, '.') - 1));

% Create new directory, change filename, and move file into new directory
mkdir(dirName);
newName = sprintf('1.%s', ext);
movefile(fin, [dirName , '/'  , newName]);

end

function copyImage(currDir, ext, dName, nCopies)
%% copyImage: copy filename n times with incrementing filename
% Go into directory fld matching filename fin
% Create n copies of fin defined by nCopies, numbered 1.ext - N.ext [ext == image extension]
cd(dName);

fin = sprintf('1.%s', ext);
for c = 2 : nCopies
    try
        fcopy = sprintf('%d.%s', c, ext);
        copyfile(fin, fcopy);
    catch e
        fprintf(2, '%s\n', e.getReport);
        cd(currDir);
    end
end

cd(currDir);
end

function crp = push2edge(im)
%% push2edge: crop image to push it's object to the right-most edge
% Determine location of object and define right-most edge of object
% Remove all columns between right-most edge of object and right edge of image
bnds     = bwboundaries(im, 'noholes');
[~, lrg] = max(cellfun(@numel, bnds));
bnd      = bnds{lrg};

% Find index at right-edge and crop
x   = bnd(:,2);
mx  = find(x == max(x), 1, 'last' );
mn  = find(x == max(x), 1 );
idx = min([x(mx) x(mn)]);
crp = im(:, 1:idx);

end

function replaceOrignalWithCrop(image_name, cropped_image)
%% replaceOrignalWithCrop: replace original image with cropped variant

try
    imwrite(cropped_image, image_name);
catch e
    fprintf(2, '%s\n', e.getReport);
end

end

function mid = phytoMorphMidline(currDir, dName)
%% phytoMorphMidline: run through customized version of phytoMorph midline algorithm
N = imageDatastore(dName);
cd(dName);

for n = 1 : numel(N.Files)
    fprintf('Running %s through custom midline algorithm...\n', N.Files{n});
    mid = sprintf('Cool data for image %d from %s', n, dName);
end

cd(currDir);
end

function midData = extractMidlineData(mid)
%% extractMidlineData: extract output from phytoMorph algorithm for a single image
% Load .mat file containing midline data
% Extract first data structure and process midline data

midData = strrep(mid, ' ', '_');
end

function resetDir(currDir, dName, ext)
%% resetDir: reset image and directory names
% Accidentally replaced backup folder with folder I screwed up
try
    cd(dName);
    fin  = sprintf('1.%s', ext);
    fout = sprintf('../%s.%s', getDirName(dName), ext);
    movefile(fin, fout);
    delete(sprintf('*.%s', ext));
    cd(currDir);
    rmdir(dName);
catch e
    fprintf(2, '%s\n', e.getReport);
    cd(currDir);
end
end


