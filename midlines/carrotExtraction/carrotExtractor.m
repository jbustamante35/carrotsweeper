function [mline, crv, smsk, pmsk, tcrd, dsts, fname] = carrotExtractor(dataIn, vis, savData, savFigs)
%% carrotExtractor: midline extraction and straightener
% This is a detailed description of this script...
%
% Note about saving with the svData parameter:
% When svData is set to true, the output is saved in a .mat file that is placed
% in the output directory [ see Usage below ]. All 4 output variables are stored
% in a single structure called CARROTS.
%
% Usage:
%   [mline, cntr, smsk, pmsk] = carrotExtractor(dataIn, vis, svData, svFigs)
%
% Input:
%   dataIn: path to directory of binary images
%   vis: boolean visualize outputs
%   svData: boolean to save output in .mat file [ see note above ]
%   svFigs: boolean to save figures of straightened carrots
%
% Output:
%   mline: cell array of midline data
%   crv: cell array of contour data
%   smsk: cell array of straightened mask images
%   pmsk: cell array of processed mask images
%   tcrd: cell array of tip coordinates
%   dsts: cell array of distance transform values along midline
%
% Usage (continued):
% If you want to save the results, the data will automatically be placed in the
% input directory as a subfolder named output_yymmdd, where 'yymmdd' corresponds
% to the year (y), month (m), and today's date (d).
%
% Example:
%   Run straightening pipeline on all images in a directory
%       dataIn  = '~/LabData/CarrotSweeper/z_datasets/masks_wi2019';
%       din     = [dataIn, '/' , 'pi-261783/binary-masks'];
%       [mline, cntr, smsk, pmsk, tcrd, dsts] = carrotExtractor(din, 1, 1, 0);
%

%% Some constants to consider playing around with
% THRESH = 300; % Minimum length to pad one or both dimensions of image
FACE = 2;   % Direction to point straightened images (original 3)

%% Load file list of binary mask images

if savData || savFigs
    dOut = sprintf('output_%s', tdate('s'));
    
    if isfolder(dataIn)
        dataOut = sprintf('%s/%s', dataIn, dOut);
    else
        dataOut = sprintf('%s/%s', fileparts(dataIn), dOut);
    end
    
    mkdir(dataOut);
end

if isfolder(dataIn)
    ext = '.png';
    img = imageDatastore(dataIn, 'FileExtensions', ext);
    
    %% Extract Midline, Contour, Straightened Image, Straightened Mask
    tot                                  = numel(img.Files);
    [mline, crv, pmsk, smsk, tcrd, dsts, fname] = deal(cell(1, tot));
    
    for n = 1 : tot
        tic;
        try
            [pmsk{n}, crv{n}, mline{n}, smsk{n}, tcrd{n}, dsts{n}] = runStraighteningPipeline(img.readimage(n));

            fname{n} = getDirName(img.Files{n});
            msg = sprintf('Successfully processed: %s', fname{n}); disp(msg);    

        catch e
            fprintf(2, 'Error in Carrot Pipeline\n%s\n', e.getReport);
        end
        toc;
    end
    
else
    img                            = imread(dataIn);
    [tot , n]                      = deal(1);
    [mline, crv, pmsk, smsk, tcrd, dsts, fname] = deal(cell(1, tot));
    
    try
        [pmsk{n}, crv{n}, mline{n}, smsk{n}, tcrd{n}, dsts{n}] = ...
            runStraighteningPipeline(img);
        
    catch e
        fprintf(2, 'Error in Carrot Extraction Pipeline\n%s\n', e.getReport);
    end
    
end

%% Show Output of processed and straightened masks
if vis
    psec = 0.7;
    figs = [];
    fnms = {};
    fnms{1} = sprintf('MidlineContourMask');
    fnms{2} = sprintf('WidthProfile');
    fnms{3} = sprintf('WidthsOnMask');
    fnms{4} = sprintf('StraightenedMask');
    
    set(figs, 'Color', 'w');
    
    for n = 1 : tot
        cla;clf;
        try
            figs = plotCarrots(n, pmsk{n}, mline{n}, crv{n}, tcrd{n}, dsts{n}, ...
                smsk{n}, psec, 0);
            
        catch e
            fprintf(2, 'Error plotting figure for data %d\n%s\n', ...
                n, e.getReport);
            
            figs = plotCarrots(n, pmsk{n}, [0 0], [0 0], [0 0], [0 0], [0 0], psec, 0);
            
        end
        
        % Save figures in output directory
        if savFigs
            [~, fName] = fileparts(dataIn);
            for fig = figs
                fnm   = sprintf('%s/%s%d_%s', dataOut, fnms{fig}, n, fName);
                savefig(fig, fnm);
                saveas(fig, fnm, 'tiffn');
            end
        end
        
    end
end

%% Save Data in output directory
% Add CSV with UID | Width | Length
if savData
    [~, fName]   = fileparts(dataIn);
    CARROTS      = v2struct(mline, crv, smsk, pmsk);
    nm           = sprintf('%s/%s_carrotExtractor_%s_%dCarrots', ...
        dataOut, tdate('s'), fName, tot);
    save(nm, '-v7.3', 'CARROTS');
    
    % Write into CSV file
    % UID [UID] | Max Width (mm) [Scale] | Length
    %% Outpt CSV with UID [UID] | Max Width (mm) [Scale] | Length
    % Extract metadata from filenames
    if isfolder(dataIn)
        % For directories of images
        tdir = dataIn;
        I    = imageDatastore(tdir);
        nms  = I.Files;
    end
    
    ID   = 'UID';
    expr = sprintf('%s_(?<id>.*?)}', ID);
    uid  = regexpi(nms, expr, 'names');
    uids = cellfun(@(x) char(x.id), uid, 'UniformOutput', 0);
    
    ID   = 'Scale';
    expr = sprintf('%s_(?<id>.*?)}', ID);
    scl  = regexpi(nms, expr, 'names');
    scls = cellfun(@(x) str2double(char(x.id)), scl, 'UniformOutput', 0);
    
    % Compute lengths and max widths
    flp_dsts       = cellfun(@(x) flipud(x), dsts, 'UniformOutput', 0)';
    flp_mln        = cellfun(@(x) flipud(x), mline, 'UniformOutput', 0)';
    [maxDst, ~]    = cellfun(@(x) max(x), flp_dsts, 'UniformOutput', 0);
    
    % Convert pix2in2mm using Scale [DPI]
    in2mm  = 25.4;
    dig    = 2;
    maxWid = cellfun(@(d,s) round((d * in2mm) / s, dig), ...
        maxDst, scls, 'UniformOutput', 0);
    maxLen = cellfun(@(w,s) round((length(w) * in2mm) / s, dig), ...
        flp_mln, scls, 'UniformOutput', 0);
    proWid = cellfun(@(x) sprintf('%.0f,', x'), dsts, 'UniformOutput', 0)';
    
    % Store UID-Wid-Len in structure format and display random carrot's data
    str  = struct('UID', uids, 'MaxWidth', maxWid, 'MaxLen', maxLen, ...
        'WidthProfile', proWid);
    
    % Convert structure to table and store as CSV
    PI   = 'PI';
    expr = sprintf('%s-(?<id>.*?)/', PI);
    pid  = regexpi(tdir, expr, 'names');
    tnm  = sprintf('%s/PI-%s.csv', dataOut, pid.id);
    tbl  = struct2table(str);
    writetable(tbl, tnm, 'FileType', 'text');
    
    tnm2 = sprintf('%s/PI-%s', dataOut, pid.id);
    writetable(tbl, tnm2, 'FileType', 'spreadsheet');
else
    % Single images [this needs to be fixed]
    tdir = dataIn;
    nms = tdir;
end

end

function figs = plotCarrots(idx, raw_mask, midline, contours, tip_crds, dsts, straight_mask, psec, f)
%% plotCarrots: plotting function for this script
% Set f to false generate figures if they don't exist
% Set f to false to overwrite existing figures
if f
    figs = 1:4;
    figs(1) = figure;
    figs(2) = figure;
    figs(3) = figure;
    figs(4) = figure;
    set(figs,  'Color',  'w');
else
    figs = 1:4;
    figs(1) = figure(1);
    figs(2) = figure(2);
    figs(3) = figure(3);
    figs(4) = figure(4);
    set(figs,  'Color',  'w');
    set(figs, 'Color', 'w');
end

%% Overlay midline, contour, tip on processed mask
fIdx = 1;
set(0, 'CurrentFigure', figs(fIdx)); fIdx = fIdx + 1;
cla;clf;

imshow(raw_mask, []);
imagesc(raw_mask);
colormap gray;
axis image;
hold on;
plt(midline, 'r-', 2);
plt(contours, 'b-', 2);
plt(tip_crds, 'g*', 5);
ttlP = sprintf('Midline and Contour on Mask\nCarrot %d', idx);
title(ttlP);

%% Bar plot of distance vector
set(0, 'CurrentFigure', figs(fIdx)); fIdx = fIdx + 1;
cla;clf;

bar(flip(dsts), 1, 'r');
ttlS = sprintf('Width Profile\nCarrot %d', idx);
title(ttlS);

%% Tick marks along midline showing widths at tick
set(0, 'CurrentFigure', figs(fIdx)); fIdx = fIdx + 1;
cla;clf;

imagesc(raw_mask);
colormap gray;
axis image;
hold on;

plt(midline, 'r.', 3);
plt(contours, 'y.', 3);
plt(tip_crds, 'g*', 6);

% Show every length / 10
lng = length(midline);
itr = ceil(lng / 15);
X   = midline(:,1) - 5;
Y   = midline(:,2) - 15;
txt = cellstr(num2str(round(dsts, 2)));

for i = 1 : itr : lng
    % Plot distance and tick marks
    text(X(i), Y(i), txt{i}, 'Color', 'b', 'FontSize', 6);
    %     text(X(i), Y(i), '+', 'Color', 'r', 'FontSize', 6); % Calibrate position
    plt(midline(i,:), 'b+', 3);
end

% Plot max distance [reverse order to read left-right]
flp_dsts       = flipud(dsts);
flp_mln        = flipud(midline);

[maxD, maxIdx] = max(flp_dsts);
maxP           = flp_mln(maxIdx,:);
maxX           = maxP(1) - 25;
maxY           = maxP(2) + 20;
maxT           = num2str(round(maxD, 2));

plt(maxP, 'k+', 8);
text(maxX, maxY, maxT, 'Color', 'k', 'FontSize', 7, 'FontWeight', 'bold');

ttlP = sprintf('Length %d pixels | Max Width %.02f pixels\nCarrot %d', ...
    lng, maxD, idx);
title(ttlP);

%% Straightened mask
set(0, 'CurrentFigure', figs(fIdx)); fIdx = fIdx + 1;
cla;clf;

flp = handleFLIP(straight_mask, 4);
imagesc(flp); colormap gray;
ttlS = sprintf('Straighted Mask\nCarrot %d', idx);
title(ttlS);

pause(psec);

end

