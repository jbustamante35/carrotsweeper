function batchStraighten(parentDir, vis, saveData, saveFigs, straightenedMasks, widthProfile, widthTicks, midlineOverlay, depth)

files = dir(parentDir);
dirFlags = [files.isdir];
subDirs = files(dirFlags);
subDirsNames = cell(1, numel(subDirs) - 2);
    for i=3:numel(subDirs)
        subDirsNames{i-2} = subDirs(i).name;
    end
% look inside of each genotype folder and run the Extractor over the binary
% masks
for i = 1:depth

    genoPath = convertStringsToChars(join([parentDir,"/",subDirsNames{i},"/"],""));
    binaryMaskPath = convertStringsToChars(join([parentDir,subDirsNames{i},"binary-masks"],"/"))

    [mline, cntr, smsk, pmsk, tcrd, dsts, fname] = carrotExtractor(binaryMaskPath, vis, saveData, saveFigs);
    disp('Midlines calculated, proceeding to write selected outputs to disk');
    if midlineOverlay == 1
        msg = sprintf('Writing midline overlays to disc (Genotype %s)', subDirsNames{i}); disp(msg);

%       %%% write overlays (how to get close cropping?)
        dirName = sprintf('%s%s', genoPath, 'midline-overlays/');
        if ~exist(dirName, 'dir')
            mkdir(dirName)
        end

        for k = 1:length(fname)
            f = figure('visible','off');
            imshow(pmsk{k}, []);
            imagesc(pmsk{k}); colormap gray; axis off;

            hold on;

            plt(mline{k}, 'r-', 2);
            plt(cntr{k}, 'b-', 2);
            plt(tcrd{k}, 'g*', 5);

            ttlP = sprintf('Midline and Contour on Mask\n%s', fname{k});
            title(ttlP);

            fnm = sprintf('%s%s%s', genoPath, 'midline-overlays/', fname{k});

            saveas(f, fnm);
            cla;clf;

        end
    end

    if straightenedMasks == 1
        msg = sprintf('Writing straightened masks to disc (Genotype %s)', subDirsNames{i}); disp(msg);

        %%% write straightened masks (they still don't look good)
        dirName = sprintf('%s%s', genoPath, 'straightened-masks/');
        if ~exist(dirName, 'dir')
            mkdir(dirName)
        end

        for k = 1:length(fname)
            straightMask = im2uint8(smsk{k});
            binarizedMask= imbinarize(straightMask);
            fnm = sprintf('%s%s%s', genoPath, 'straightened-masks/', fname{k});
            imwrite(rot90(binarizedMask), fnm, 'PNG');
        end
    end


    if widthProfile == 1
        msg = sprintf('Writing width profiles to disc as .csv files (Genotype %s)', subDirsNames{i}); disp(msg);

        %%% write width profiles (requires r2019a)
        dirName = sprintf('%s%s', genoPath, 'width-vectors/');
        if ~exist(dirName, 'dir')
            mkdir(dirName)
        end

        for k = 1:length(fname)
            [filepath,name,ext] = fileparts(fname{k});
            fnm = sprintf('%s%s%s%s', genoPath, 'width-vectors/', name, '.csv');
            writematrix(dsts{k}, fnm);
        end
    end


	if widthProfile == 2
        sprintf('Writing width profiles to disc as barplots (Genotype %s)', subDirsNames{i})


        dirName = sprintf('%s%s', genoPath, 'width-profiles/');
        if ~exist(dirName, 'dir')
            mkdir(dirName)
        end

        for k = 1:length(fname)
            f = figure('visible','off');

            bar(flip(dsts{k}), 1, 'r');
            ttlS = sprintf('Width Profile\nRoot %d', k);
            title(ttlS);

            fnm = sprintf('%s%s%s', genoPath, 'width-profiles/', fname{k});
            saveas(f, fnm);
            cla;clf;
        end
    end

	if widthTicks == 1
        dirName = sprintf('%s%s', genoPath, 'width-slices/');
        if ~exist(dirName, 'dir')
            mkdir(dirName)
        end

        for k = 1:length(fname)
            % Still need to do this
        end
    end


end

end

