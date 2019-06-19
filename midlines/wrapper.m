%%
% Get list of genotype folders to process

parentDir = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/images/model-validation/2019-groundtruth/original';

files = dir(parentDir);
dirFlags = [files.isdir];
subDirs = files(dirFlags);
subDirsNames = cell(1, numel(subDirs) - 2);
    for i=3:numel(subDirs)
        subDirsNames{i-2} = subDirs(i).name;
    end
    
%%

% Set variables for extractor
vis = 0;
saveData = 0;
saveFigs = 0;

% Set variables for wrapper
straightenedMasks = 0; % smsk
widthProfile = 'barplot'; % dsts = 'vector' or barplot = 'barplot'
widthTicks = 0; % normals
midlineOverlay = 1; % pmsk with mline and tcrd
    
%%

% look inside of each genotype folder and run the Extractor over the binary
% masks
for i = 1:length(subDirsNames)
    
    genoPath = convertStringsToChars(join([parentDir,"/",subDirsNames{i},"/"],""))
    binaryMaskPath = convertStringsToChars(join([parentDir,subDirsNames{i},"binary-masks"],"/"))
    [mline, cntr, smsk, pmsk, tcrd, dsts, fname] = carrotExtractor(binaryMaskPath, vis, saveData, saveFigs);

    if midlineOverlay == 1
%       %%% write overlays (how to get close cropping?)
        cd(genoPath);
        mkdir('midline-overlays');
        
        for k = 1:length(fname)   
            imshow(pmsk{k}, []);
            imagesc(pmsk{k}); colormap gray; axis off;

            hold on;

            plt(mline{k}, 'r-', 2);
            plt(cntr{k}, 'b-', 2);
            plt(tcrd{k}, 'g*', 5);
            ttlP = sprintf('Midline and Contour on Mask\nRoot %d', k);
            title(ttlP);

            fnm   = sprintf('%s%s%s', genoPath, 'midline-overlays/', fname{k});
            saveas(figure(1), fnm);
            cla;clf;

        end
    end
    
    if straightenedMasks == 1
        %%% write straightened masks (they still don't look good)
        cd(genoPath);
        mkdir('straightened-masks');
        cd 'straightened-masks'

        for k = 1:length(fname)   
            straightMask = im2uint8(smsk{k});
            binarizedMask= imbinarize(straightMask);
            imwrite(rot90(binarizedMask), fname{k}, 'PNG');
        end
    end
    
    
    if widthProfile == 'vector'
        %%% write width profiles (requires r2019a)
        cd(genoPath);
        mkdir('width-profiles');
        cd 'width-profiles'

        for k = 1:length(fname)   
            writematrix(dsts{k}, fname{k});
        end
    end
    
    
	if widthProfile == 'barplot'
        %%% write width profiles (requires r2019a)

    end 
end

%%
        
