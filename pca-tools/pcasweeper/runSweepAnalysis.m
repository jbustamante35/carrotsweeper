function D = runSweepAnalysis(BW, numX, numY, stps, vis, sav)
%% runSweepAnalysis: extract contours from dataset, run PCA, then sweep through PCs
% This function runs a pipeline in which the user loads a cell array of bw
% images, then defines the number of Principal Components (PCs) to extract
% from PCA on the x-/y-coordinates. The number of respective PCs to use for
% the x and y coordinate are set by the numX and numY parameters.
%
% The data from PCA then go through the sweeping function, where each PC from
% the x-/y-coordinates are iteratively swept up or down by a standard
% deviation, up to the number of steps defined by the stps parameter.
%
% Additional data can be visualized if the boolean vis parameter is set to 1;
% otherwise only the figures from the PC sweep will be shown. The boolean sv
% parameter saves the figure handles as both .fig and .tif files, and the PCA
% output and PC sweep output as a single .mat file.
%
% Usage:
%   D = runSweepAnalysis(BW, numX, numY, stps, vis, sav)
%
% Input:
%   BW: cell array of bw images to be analyzed
%   numX: number of PCs to extract for PCA on the x-coordinates
%   numY: number of PCs to extract for PCA on the y-coordinates
%   stps: number of steps above/below standard deviations to sweep through
%   vis: boolean to show only sweep figures (0) or additional figures (1)
%   sav: boolean to save figures as .fig and .tif files and data as a .mat file
%
% Output:
%   D: structure containing PCA data and PC sweep data
%
% Examples:
%   dat = '/mnt/tetra/JulianBustamante/HypoQuantyl/Contours';
%   C   = load(sprintf('%s/180830_scott/180830_carrotPCA.mat', dat)); % 700 roots unaligned dataset
%   C   = C.CARROT;
%   BW  = C.bw';
%   D   = runSweepAnalysis(BW, 5, 5, 7, 0, 0)
%
%   dat = '/mnt/tetra/JulianBustamante/HypoQuantyl/Contours';
%   V   = load(sprintf('%s/180830_scott/carrotPCAb.mat', dat)); % 328 roots aligned dataset
%   BW  = V.bw';
%   W   = runSweepAnalysis(BW, 3, 2, 5, 1, 1)
%

%% Prepare empty figure plots
if vis
%     figs = [];
%     fnms = {};

%     figs(1) = figure; % Reconvert simulated contours
%     figs(2) = figure; % Mean contour vs random contours
%     figs(3) = figure; % All Sweeps on Mean Contour
%     figs(4) = 4; % Sweep through x-coordinate PCs
%     figs(5) = 5; % Sweep through y-coordinate PCs
    [figs , fnms] = makeBlankFigures(5, 1);

    fnms{1} = sprintf('%s_ReconvertedOnImage', tdate);
    fnms{2} = sprintf('%s_MeanVsRandomContours', tdate);
    fnms{3} = sprintf('%s_SweepsOnMeanContour', tdate);
    fnms{4} = sprintf('%s_PCSweep_xCoords', tdate);
    fnms{5} = sprintf('%s_PCSweep_yCoords', tdate);

else
    figs = 1 : 5; % PC sweep figures are figs(4:5)

    fnms{1} = '';
    fnms{2} = '';
    fnms{3} = '';
    fnms{4} = sprintf('%s_PCSweep_xCoords', datestr(now, 'yymmdd'));
    fnms{5} = sprintf('%s_PCSweep_yCoords', datestr(now, 'yymmdd'));

end

%% Extract and Rasterize normalized contours from BW images
[X, Y, CNTR] = extractAndRasterize(BW);

%% Run PCA on x-/y-coordinates
px = pcaAnalysis(X, numX, 0, 'xCoords');
py = pcaAnalysis(Y, numY, 0, 'yCoords');

%% Sweep through all PCs
[scFull, smFull] = ...
    performSweep('pcaX', px, 'pcaY', py, 'nsteps', stps, 'figs', 1:2, 'sav', 0);

%% Set x-/y-limits equal [figure out how to set this dynamically]
% xl = [-900 200];
% yl = [-150 150];
% xl = getMax(px);
% yl = getMax(py);

% x-coordinates
% figclr(4);
% for s = 1 : numX
%     subplot(2, 2, s);
%     %     xlim(xl);
%     %     ylim(yl);
% end
%
% y-coordinates
% figclr(5);
% for s = 1 : numY
%     subplot(2, 1, s);
%     %     xlim(xl);
%     %     ylim(yl);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot re-converted simulated contour onto bw image
if vis
    figclr(3);

    numCV = 16;
    rows  = 4;
    cols  = 4;
    for p = 1 : numCV
        % Get random index and extract data
        rIdx = pullRandom(BW);
        img  = BW{rIdx};
%         apt  = CNTR{rIdx}.getAnchorPoint;
%         aid  = CNTR{rIdx}.getAnchorIndex;

        % Convert normalized input and simulated contour to raw image coordinates
        cInp = [px.InputData(rIdx,:) ; py.InputData(rIdx,:)]';
%         cInp = norm2raw(nInp, apt, aid);
        cSim = [px.SimData(rIdx) ; py.SimData(rIdx)]';
%         cSim = norm2raw(nSim, apt, aid);

        % Overlay converted input and simulated contour on bw image
        subplot(rows, cols, p);        
        imagesc(img);
        hold on;
        plt(cInp, 'g-', 1);
        plt(cSim, 'y-', 1);

        colormap gray;
        axis ij;
        axis tight;
        ttl = sprintf('Contour %d', rIdx);
        title(ttl);
    end

    %% Compare mean contour with actual contours
    figclr(4);    
    
    CT    = cat(1, CNTR{:});
    numCT = ceil(numel(CNTR) / 10);
    for c = 1 : numCT
        rcntr = CT(pullRandom(CT)).NormalizedOutline;
        plt(rcntr, 'm--', 1);
        hold on;
    end

    mean_contour = smFull{1}{1}.mean;
    plt(mean_contour, 'k-', 5);

    axis ij;
    ttl = sprintf('%d contours on mean contour', numCT);
    title(ttl);

    %% Plot all PC sweeps on single mean contour
    figclr(5);
    hold on;

    for d = 1 : size(smFull, 1)
        dim = smFull(d,:);
        for p = 1 : size(dim, 2)
            pc = dim{p};
            for s = 1 : size(pc, 2)
                try
                    swp = pc{s};
                    plt(swp.up, 'g--', 1);
                    plt(swp.down, 'r--', 1);
                catch
                    continue;
                end
            end
        end
    end

    plt(mean_contour, 'k-', 5);

    ttl = sprintf('%d Swept PCs on mean contour', stps);
    title(ttl);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Figures and Dataset
D  = v2struct(px, py, scFull, smFull, CNTR);

if sav
    saveFiguresJB(figs, fnms, 0);

    nm = sprintf('%s_carrotPCA_analysis_%dCarrots', tdate, numel(CNTR));
    save(nm, '-v7.3', 'D');
end
end


function [X, Y, C] = extractAndRasterize(BW)
%% extractAndRasterize: create ContourJB then split and rasterize coordinates
sz         = 800;
[CNTR , C] = cellfun(@(x) extractContour(x, sz, 'Normalized'), ...
    BW, 'UniformOutput', 0);

bndX = cellfun(@(x) getDim(x,1), CNTR, 'UniformOutput', 0);
bndY = cellfun(@(x) getDim(x,2), CNTR, 'UniformOutput', 0);
X    = rasterizeImagesHQ(bndX);
Y    = rasterizeImagesHQ(bndY);

end

