function [mline, crv, smsk, pmsk, tcrd, dsts, fname] = carrotStraightener(dataIn, savData, savFigs, vis, par)
%% carrotStraightener: midline extraction and straightener
% This is a detailed description of this script...
%
% Note about saving with the svData parameter:
% When svData is set to true, the output is saved in a .mat file that is placed
% in the output directory [ see Usage below ]. All 4 output variables are stored
% in a single structure called CARROTS.
%
% Usage:
%   [mline, crv, smsk, pmsk, tcrd, dsts, fname] = ...
%         carrotStraightener(dataIn, savData, savFigs, vis, par)
%
% Input:
%   dataIn: path to directory of binary images
%   vis: boolean visualize outputs
%   svData: boolean to save output in .mat file [ see note above ]
%   svFigs: boolean to save figures of straightened carrots
%   par: boolean to run on single thread (0) or with parallel processing (1)
%
% Output:
%   mline: cell array of midline data
%   crv: cell array of contour data
%   smsk: cell array of straightened mask images
%   pmsk: cell array of processed mask images
%   tcrd: cell array of tip coordinates
%   dsts: cell array of distance transform values along midline
%   fname: cell array of filenames of images
%
% Usage (continued):
% If you want to save the resulzts, the data will automatically be placed in the
% input directory as a subfolder named output_yymmdd, where 'yymmdd' corresponds
% to the year (y), month (m), and today's date (d).
%
% Example:
%   Run pipeline on all images in a directory with paraellelization
%       dataIn  = '~/LabData/CarrotSweeper/z_datasets/masks_wi2019';
%       din     = [dataIn, '/' , 'pi-261783/binary-masks'];
%       [mline, cntr, smsk, pmsk, tcrd, dsts, fname] = ...
%           carrotStraightener(din, 1, 1, 1, 1);
%

%% Load file list of binary mask images
% Save outputted data in respective directories
if savData || savFigs
    dMat = sprintf('output-%s', tdate('s'));
    dOvr = 'mask-overlays';
    dWid = 'width-profiles';
    dNrm = 'normal-overlays';
    dMsk = 'straight-masks';
    
    if isfolder(dataIn)
        matOut = sprintf('%s/%s', fileparts(dataIn), dMat);
        ovrOut = sprintf('%s/%s', fileparts(dataIn), dOvr);
        widOut = sprintf('%s/%s', fileparts(dataIn), dWid);
        nrmOut = sprintf('%s/%s', fileparts(dataIn), dNrm);
        mskOut = sprintf('%s/%s', fileparts(dataIn), dMsk);
    else
        matOut = sprintf('%s/%s', fileparts(fileparts(dataIn)), dMat);
        ovrOut = sprintf('%s/%s', fileparts(fileparts(dataIn)), dOvr);
        widOut = sprintf('%s/%s', fileparts(fileparts(dataIn)), dWid);
        nrmOut = sprintf('%s/%s', fileparts(fileparts(dataIn)), dNrm);
        mskOut = sprintf('%s/%s', fileparts(fileparts(dataIn)), dMsk);
    end
    
    mkdir(matOut);
    mkdir(ovrOut);
    mkdir(widOut);
    mkdir(nrmOut);
    mkdir(mskOut);
    
    % Open figures if none are open
    nf   = 4;
    figs = 1 : nf;
    if isempty(findobj('type', 'figure')) || ...
            numel(findobj('type', 'figure')) < nf
        % Generate n figures        
        for n = 1 : nf
            figs(n) = figure;
        end
        set(figs, 'Color', 'w');
    end
    
end

if isfolder(dataIn)
    ext = '.png';
    img = imageDatastore(dataIn, 'FileExtensions', ext);
    
    %% Extract Midline, Contour, Straightened Image, Straightened Mask
    tot                                               = numel(img.Files);
    [mline, crv, pmsk, smsk, tcrd, dsts, fname, nrms] = deal(cell(1, tot));
    
    if par
        %% Run through images with parallel processing
        tt = tic;
        parfor n = 1 : tot
            t = tic;
            try
                fname{n} = getDirName(img.Files{n});
                fprintf('\n%s\n', repmat('=', 1, 80));
                fprintf('Processing image %d of %d\n%s', n, tot, fname{n});
                fprintf('\n%s\n', repmat('-', 1, 80));
                
                [pmsk{n}, crv{n}, mline{n}, smsk{n}, tcrd{n}, dsts{n}, nrms{n}] = ...
                    runStraighteningPipeline(img.readimage(n));
                
                fprintf('\n%s\n', repmat('-', 1, 80));
                fprintf('Successfully processed %s\n', fname{n});
                
            catch e
                fprintf(2, 'Error processing %s\n%s\n', fname{n}, e.getReport);
            end
            fprintf('Pipeline finished in %.02f sec', toc(t));
            fprintf('\n%s\n', repmat('-', 1, 80));
        end
        fprintf('DONE! [%.02f sec]\n', toc(tt));
    else
        %% Run through images with normal for loop
        tt = tic;
        for n = 1 : tot
            t = tic;
            try
                fname{n} = getDirName(img.Files{n});
                fprintf('\n%s\n', repmat('=', 1, 80));
                fprintf('Processing image %d of %d\n%s', n, tot, fname{n});
                fprintf('\n%s\n', repmat('-', 1, 80));
                
                [pmsk{n}, crv{n}, mline{n}, smsk{n}, tcrd{n}, dsts{n}, nrms{n}] = ...
                    runStraighteningPipeline(img.readimage(n));
                
                fprintf('\n%s\n', repmat('-', 1, 80));
                fprintf('Successfully processed %s\n', fname{n});
                
            catch e
                fprintf(2, 'Error processing %s\n%s\n', fname{n}, e.getReport);
            end
            fprintf('Pipeline finished in %.02f sec', toc(t));
            fprintf('\n%s\n', repmat('=', 1, 80));
        end
        fprintf('DONE! [%.02f sec]\n', toc(tt));
    end
else
    %% Run on single image
    t = tic;
    fname{1}                                   = getDirName(dataIn);
    img                                        = imread(dataIn);
    [tot , n]                                  = deal(1);
    [pmsk, crv, mline, smsk, tcrd, nrms, dsts] = deal(cell(1, tot));
    
    try        
        fprintf('%s\n', repmat('=', 1, 80));
        fprintf('Processing image %d of %d\n%s', n, tot, fname{n});        
        fprintf('%s\n', repmat('-', 1, 80));
        [pmsk{n}, crv{n}, mline{n}, smsk{n}, tcrd{n}, dsts{n}, nrms{n}] = ...
            runStraighteningPipeline(img);
        fprintf('Pipeline finished in %.02f sec', toc(t));
        fprintf('%s\n', repmat('=', 1, 80));
    catch e
        fprintf(2, 'Error in Carrot Extraction Pipeline\n%s\n', e.getReport);
    end
    
end

%% Show Output of processed and straightened masks
if vis
    psec = 0.7;
    car  = 'carrots'; % For making cleaner figure titles
    
    for n = 1 : tot
        nm    = fixtitle(fname{n}, car);
        nmidx = sprintf('%03d of %03d', n, tot);
        fnm   = {nm , nmidx};
        fIdxs = plotCarrots(fnm, pmsk{n}, mline{n}, crv{n}, tcrd{n}, ...
            dsts{n}, nrms{n}, smsk{n}, psec, figs);
        
        %% Save figure in it's own directory named as the filename
        if savFigs
            % Contour-Midline-Mask overlay
            ovrFig = sprintf('%s/%s', ovrOut, fname{n});
            saveas(fIdxs(1), ovrFig, 'png');
            
            % Width Profile
            widFig = sprintf('%s/%s', widOut, fname{n});
            saveas(fIdxs(2), widFig, 'png');
            
            % Widths and Normals along Midline
            nrmFig = sprintf('%s/%s', nrmOut, fname{n});
            saveas(fIdxs(3), nrmFig, 'png');
        end
    end
end

%% Save Data in output directory
if savData
    % Save .mat and .csv files in output-$date$ directory
    if isfolder(dataIn)
        [~, fName] = fileparts(fileparts(dataIn));
    else
        [~, fName] = fileparts(dataIn);
    end
    
    flds    = {'fieldNames', 'mline', 'crv', 'smsk', 'pmsk', ...
        'tcrd', 'dsts', 'fname'};
    CARROTS = v2struct(flds);
    nm      = sprintf('%s/%s_carrotExtractor_%s_%dCarrots', ...
        matOut, tdate('s'), fName, tot);
    save(nm, '-v7.3', 'CARROTS');
    
    % Save straightened masks in straight-masks directory
    cellfun(@(im,nm) imwrite(im, [mskOut '/' nm], 'png'), ...
        smsk, fname, 'UniformOutput', 0);
    
    %% Outpt CSV with UID [UID] | Max Width (mm) [Scale] | Length
    % Extract metadata from filenames
    if isfolder(dataIn)
        % For directories of images
        tdir = dataIn;
        I    = imageDatastore(tdir);
        nms  = I.Files;
    else
        nms = dataIn;
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
    flp_dsts    = cellfun(@(x) flipud(x), dsts, 'UniformOutput', 0)';
    flp_mln     = cellfun(@(x) flipud(x), mline, 'UniformOutput', 0)';
    [maxDst, ~] = cellfun(@(x) max(x), flp_dsts, 'UniformOutput', 0);
    
    % Convert pix2in2mm using Scale [DPI]
    in2mm  = 25.4; % Convert inches to millimeters
    dig    = 2;    % Round to n digits
    maxWid = cellfun(@(d,s) round((d * in2mm) / s, dig), ...
        maxDst, scls, 'UniformOutput', 0);
    maxLen = cellfun(@(w,s) round((length(w) * in2mm) / s, dig), ...
        flp_mln, scls, 'UniformOutput', 0);
    proWid = cellfun(@(x) sprintf('%.0f,', x'), dsts, 'UniformOutput', 0)';
    
    % Store UID-Wid-Len in structure format and display random carrot's data
    str  = struct('UID', uids, 'MaxWidth', maxWid, 'MaxLen', maxLen, ...
        'WidthProfile', proWid);
    
    % Convert structure to table and store as CSV
    if isfolder(dataIn)
        [~, idDir] = fileparts(fileparts(dataIn));
    else
        [mskPath, ~] = fileparts(fileparts(dataIn));
        [idPath, ~]  = fileparts(mskPath);
        [~, idDir]   = fileparts(idPath);
    end
    
    tnm1 = sprintf('%s/%s.csv', matOut, idDir);
    
    tbl  = struct2table(str);
    writetable(tbl, tnm1, 'FileType', 'text');
    
    tnm2 = sprintf('%s/%s', matOut, idDir);
    writetable(tbl, tnm2, 'FileType', 'spreadsheet');
    
else
    % Single images [this needs to be fixed]
    tdir = dataIn;
    nms  = tdir;
end

end

function figs = plotCarrots(fnmidx, raw_mask, midline, contours, tip_crds, dsts, nrms, straight_mask, psec, figs)
%% plotCarrots: plotting function for this script
%
%
%

%% Overlay midline, contour, tip on processed mask
fIdx  = 1;
fname = fnmidx{1};
fnidx = fnmidx{2};

try
    figclr(fIdx);    
    myimagesc(raw_mask, 'gray', 'image', 'on');
    hold on;
    plt(midline, 'r-', 2);
    plt(contours, 'b-', 2);
    plt(tip_crds, 'g*', 5);
    ttlP = sprintf('Midline and Contour on Mask [%s]\n%s', fnidx, fname);
    title(ttlP, 'FontSize', 10);
    
catch
    figname = 'Mask Overlay';
    fprintf(2, 'Error plotting figure %d for %s [%s]\n', fIdx, fname, figname);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bar plot of distance vector
fIdx = fIdx + 1;

try
    figclr(fIdx);
    
    bar(flip(dsts), 1, 'r');
    axis image;
    ttlS = sprintf('Width Profile [%s]\n%s', fnidx, fname);
    title(ttlS, 'FontSize', 10);
    
catch
    figname = 'Width Profile';
    fprintf(2, 'Error plotting figure %d for %s [%s]\n', fIdx, fname, figname);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tick marks along midline showing widths and normals at tick
fIdx = fIdx + 1;

try
    figclr(fIdx);        
    myimagesc(raw_mask, 'gray', 'image', 'on');
    hold on;
    
    plt(midline, 'r-', 2);
    plt(contours, 'y-', 2);
    plt(tip_crds, 'g*', 5);
    
    % Show every length / 10
    lng = length(midline);
    itr = ceil(lng / 15);
    xos = 5;  % x-offset
    yos = 15; % y-offset
    X   = midline(:,1) - xos;
    Y   = midline(:,2) - yos;
    txt = cellstr(num2str(round(dsts', 2)));
    
    % Extract indices for normal vectors
    mIdxs = 1 : itr : lng;
    dscl  = ceil(size(raw_mask,1) / 2) + 1;
    
    nrmO = nrms.OuterData.eCrds;
    nrmI = nrms.InnerData.eCrds;
    for i = mIdxs
        % Plot distance and tick marks
        text(X(i), Y(i), txt{i}, 'Color', 'b', 'FontSize', 6);
        %     text(X(i), Y(i), '+', 'Color', 'r', 'FontSize', 6); % Calibrate position
        nIdx = extractIndices(i, dscl);
        plt(nrmO(nIdx,:), 'r-', 1);
        plt(nrmI(nIdx,:), 'b-', 1);
        plt(midline(i,:), 'b+', 3);
    end
    
    % Plot max distance
    [maxD, maxIdx] = max(dsts);
    maxP           = midline(maxIdx,:);
    maxX           = maxP(1) - 25;
    maxY           = maxP(2) + 20;
    maxT           = num2str(round(maxD, 2));
    
    % Figure title
    ttlP = sprintf('Length %d pixels | Max Width %.0f pixels [%s]\n%s', ...
        lng, maxD, fnidx, fname);
    title(ttlP, 'FontSize', 10);
    
    % Plot max normal too
    % Issues when max width is right at base [not enough envelope coordinates]
    mIdx = extractIndices(maxIdx, dscl);
    plt(nrmO(mIdx,:), 'm-', 1);
    plt(nrmI(mIdx,:), 'g-', 1);
    
    plt(maxP, 'k+', 8);
    text(maxX, maxY, maxT, 'Color', 'k', 'FontSize', 7, 'FontWeight', 'bold');
    
catch
    figname = 'Normals on Mask';
    fprintf(2, 'Error plotting figure %d for %s [%s]\n', fIdx, fname, figname);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Straightened mask
fIdx = fIdx + 1;

try
    figclr(fIdx);
    
    flp  = handleFLIP(straight_mask, 3);
    myimagesc(flp, 'gray', 'image', 'on');
    
    ttlS = sprintf('Straighted Mask [%s]\n%s', fnidx, fname);
    title(ttlS, 'FontSize', 10);
    
catch
    figname = 'Straightened Mask';
    fprintf(2, 'Error plotting figure %d for %s [%s]\n', fIdx, fname, figname);
    
end

pause(psec);

end

