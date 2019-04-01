function [out, tipBest, tipAll] = isolate_carrot_Roots(msk, vis, DEBUG)
%% isolate_carrot_Roots: customized version of PhytoMorph's isolate_Roots
% This function is difficult to understand...
%
% Usage:
%   [out, tipAll] = isolate_carrot_Roots(msk, vis, DEBUG)
%
% Input:
%   msk: binary mask image
%   vis: boolean to visualize output
%   DEBUG: boolean to run debug mode to find best values for constants
%
% Output:
%   out: output structure containing midline and contour data
%   tipAll: structure containing data from iterating through range of constants
%

%% Constants and Data structure set-up
% Contstants for finding highest curvature
KSNIP   = 163; % Original 50
SMOOTH1 = 5;   % Original 15
SMOOTH2 = 10;  % Original 30

% Thresholds for filtering out small sizes of data
CURV_THRESH = 80;  % Original 80
SKEL_THRESH = 300; % Original 500
MASK_THRESH = 300; % Original 300

% Constants for tracking point at gradient
MAX_STEP = 50000;     % Original 50000
RHO      = 20;       % Original 20
RAD      = pi / 2;   % Original pi / 2
PDENSITY = [20 200]; % Original [20 200]
WSIGMA   = 0.3;      % Original 0.3

% Output data structures
crv      = struct('level', [], 'data', [], 'length', []);
midlines = struct('data', []);
out      = struct('midlines', [], 'contours', []);

%% Pre-process mask by facing object left-right
try
    if size(msk, 3) > 1
        msk = rgb2gray(msk);
    end

    bak        = msk;
    msk(:,end) = [];
    msk        = handleFLIP(msk, []);
    msk        = padarray(msk, [0 MASK_THRESH], 'pre', 'replicate');
    %     Io         = msk;
    %
    %     %% Get skeleton structure and store contour
    %     % This is where everything goes wrong [ often creates empty data ]
    %     % Perhaps playing around with the SKEL_THRESH value will help?
    %     msk      = Io;
    %     gI       = double(imcomplement(bwdist(msk)));
    %     gI       = imfilter(gI,fspecial('gaussian', [31 31], 7));
    %     [g1, g2] = gradient(gI);
    %     thresh   = graythresh(msk);
    %     skel     = msk < thresh;
    %     skel     = bwareaopen(skel, SKEL_THRESH);
    %
    %     % Setting edge with object to 0
    %     if sum(~msk(:, end)) > 0
    %         % Right side
    %         cIdx = size(msk, 2);
    %     else
    %         % Left sie
    %         cIdx = 1;
    %     end
    %
    %     %
    %     saveVec       = skel(:, cIdx);
    %     skel(:, cIdx) = 0;
    %
    %     %
    %     cB          = imclearborder(skel);
    %     cB(:, cIdx) = saveVec;
    %     skel        = cB;
    %     skel        = imfill(skel, 'holes');
    %
    %     %
    %     cB   = imclearborder(skel);
    %     skel = skel - cB;
    %     skel = imclose(skel, strel('disk', 21, 0));
    %     skel = imfill(skel, 'holes');
    %
    %     %
    %     C = contourc(double(skel), [1 1]);
    %
    %     %% Get the crv structure
    %     str     = 1;
    %     c       = 1;
    %     ttlCrds = size(C, 2);
    %
    %     while str < ttlCrds
    %         ed            = str + C(2, str);
    %         crv(c).level  = C(1, str);
    %         crv(c).data   = C(:, str+1:ed);
    %         crv(c).length = size(crv(c).data, 2);
    %         c             = c + 1;
    %         str           = ed + 1;
    %     end
    %
    %     % Remove curves below threshold
    %     crv(cell2mat(arrayfun(@(x) x.length < CURV_THRESH, ...
    %         crv, 'UniformOutput', 0))) = [];

    %% Find the highest curvature
    if DEBUG
%         %% Debug curvature finder with differing values of constants
%         % A main point of failure is when the function to compute the point of
%         % highest curvature identifies a coordinate on the side of the object,
%         % rather than the desired tip of the structure. This debug function iterates
%         % through a broad range of values for the following constants and plots
%         % the highest curvature point.
% 
%         %% Clear figure and show mask
%         cla;clf;
%         imshow(msk, []);
%         hold on;
% 
%         %% Set debug state
%         minItr  = 1;
%         maxItr  = 151;
%         stpItr  = 10;
%         rng     = minItr : round(maxItr / stpItr) : maxItr;
%         halfIdx = ceil(size(crv.data, 2) / 2);
%         halfCrd = crv.data(:, halfIdx);
%         f       = 0;
% 
%         % Dewitt, Anakin
%         [tipIdx, tipCrds, tipBest, tipAll] = runDebugMode(crv, rng, halfCrd, f);
% 
%         %%
%         img     = double(bwdist(~skel));
%         g1      = -g1;
%         g2      = -g2;
% 
%         %
%         t1      = ba_interp2(g1, tipCrds(1), tipCrds(2));
%         t2      = ba_interp2(g2, tipCrds(1), tipCrds(2));
% 
%         %
%         N     = [t2 t1];
%         N     = -N / norm(N);
%         T     = [N(2) -N(1)];
%         iDirc = [T ; N];
% 
%         %%
%         midlines.data = trackFromPointAtGradient_carrot(img, tipCrds, iDirc, ...
%             MAX_STEP, RHO, RAD, PDENSITY, WSIGMA);
% 
    else
        %% Use my method for extracting contours
        CNTR_LENGTH = 4000;
        skel        = ~msk;
        %         msk(:,1:10) = 1;
        cntr = extractContour(skel, CNTR_LENGTH);
        crv  = struct('level', size(cntr.InterpOutline', 1), ...
            'data', cntr.InterpOutline', ...
            'length', size(cntr.InterpOutline', 2));

        %% Run with default constant values
        try
            oInit        = cwtK(crv.data', SMOOTH1);
            [~, tipIdx]  = min(oInit.K);
            oFine        = cwtK(crv.data', SMOOTH2);
            [~, fineIdx] = min((oFine.K(tipIdx - KSNIP : tipIdx + KSNIP)));
            tipIdx       = tipIdx + (fineIdx - KSNIP - 1);

        catch e
            fprintf(2, 'Reverting to original constants...\n%s\n', e.message);
            KSNIP   = 50;
            SMOOTH1 = 15;
            SMOOTH2 = 30;

            % Try again with original constant values
            oInit        = cwtK(crv.data', SMOOTH1);
            [~, tipIdx]  = min(oInit.K);
            oFine        = cwtK(crv.data', SMOOTH2);
            [~, fineIdx] = min((oFine.K(tipIdx - KSNIP : tipIdx + KSNIP)));
            tipIdx       = tipIdx + (fineIdx - KSNIP - 1);

        end

        %
        gI       = double(bwdist(msk));
        gI       = imfilter(gI,fspecial('gaussian', [31 31], 7));
        [g1, g2] = gradient(gI);

        %
        img = double(bwdist(msk));
        g1  = -g1;
        g2  = -g2;

        %
        tipCrds = crv.data(:, tipIdx);
        t1      = ba_interp2(g1, tipCrds(1), tipCrds(2));
        t2      = ba_interp2(g2, tipCrds(1), tipCrds(2));

        %
        N     = [t2 t1];
        N     = -N / norm(N);
        T     = [N(2) -N(1)];
        iDirc = [T ; N];

        %% Trace midline from tip coordinate to
        [tipBest, tipAll] = deal([]);
        midlines.data     = trackFromPointAtGradient_carrot(img, tipCrds, iDirc, ...
            MAX_STEP, RHO, RAD, PDENSITY, WSIGMA);

    end

    %% Visualize outputs
    if vis
        cla;clf;
        imshow(Io, []);
        hold on;

        plt(crv.data', 'r.', 5);
        plt(midlines.data', 'g--', 2);
        plt(crv.data(:, tipIdx)', 'g*', 8);

    end

    %% Output final results
    out.midlines = midlines;
    out.contours = crv;

catch e
    fprintf(2, 'Error isolating root\n%s\n', e.getReport);
    out.ME = e;
end

end

function [tipIdx, tipCrds, tipBest, tipAll] = runDebugMode(crv, rng, halfCrd, f)
%% Make color array and set-up output structure
% Iterate through colors each time base for loop iterates
cAll = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'];
cIdx = [ 1    2    3    4    5         7      ];
cols = cAll(cIdx);
cItr = 1;

% Cutoff percentage to stop loops when near ground truth point
X        = zeros(1,1);
tipAll   = struct('KSNIP', X, 'SMOOTH1', X, 'SMOOTH2', X, ...
    'idx', X, 'crds', [], 'pct', X);
tipBest  = struct('KSNIP', 0, 'SMOOTH1', 0, 'SMOOTH2', 0, ...
    'idx', 0, 'crds', [], 'pct', []);
bestPct = 1.00;
allN  = 1;

%% Iterate through range of constant values
% TODO Fix this triple `for` loop I hate it
for k = rng
    for kk = rng
        for kkk = rng
            % Set constant values and set up parameters for output
            KSNIP         = k;   % Original 50
            SMOOTH_VALUE  = kk;  % Original 15
            SMOOTH_VALUE2 = kkk; % Original 30

            col = [cols(cItr) 'x'];
            err = 1;

            % Run algorithm to find point of highest curvature
            try
                oInit        = cwtK(crv.data', {SMOOTH_VALUE});
                [~, tipIdx]  = min(oInit.K);
                oFine        = cwtK(crv.data', {SMOOTH_VALUE2});
                [~, fineIdx] = min((oFine.K(tipIdx - KSNIP : tipIdx + KSNIP)));
                tipIdx       = tipIdx + (fineIdx - KSNIP - 1);
                tipCrds      = crv.data(:, tipIdx);
            catch
                tipCrds = crv.data(:,1);
                col     = 'yo';
                err     = 2;
            end

            % Store values closest to ground-truth
            currPct = abs(1 - (tipCrds \ halfCrd));
            if currPct <= bestPct
                bestPct         = currPct;
                tipBest.KSNIP   = KSNIP;
                tipBest.SMOOTH1 = SMOOTH_VALUE;
                tipBest.SMOOTH2 = SMOOTH_VALUE2;
                tipBest.idx     = tipIdx;
                tipBest.crds    = tipCrds;
                tipBest.pct     = currPct;
            end

            % Overlay tip coordinate over mask
            if f ~= 0
                set(0, 'CurrentFigure', f);
                plt(tipCrds', col, 8);
                drawnow;
            else
                tipAll.KSNIP(allN)   = KSNIP;
                tipAll.SMOOTH1(allN) = SMOOTH_VALUE;
                tipAll.SMOOTH2(allN) = SMOOTH_VALUE2;
                tipAll.idx(allN)     = tipIdx{1};
                tipAll.crds(:, allN) = tipCrds;
                tipAll.pct(allN)     = currPct;
                allN                 = allN + 1;
            end

            % Output constant values and percentages away from true tip
            fprintf(err, ...
                'CURR %.02f | BEST %.02f | K %d | S1 %d | S2 %d\n', ...
                currPct, bestPct, KSNIP, SMOOTH_VALUE, SMOOTH_VALUE2);
        end
    end

    % Cycle through colors after base for-loop iterates
    if cItr >= numel(cols)
        cItr = 1;
    else
        cItr = cItr + 1;
    end
end

%% Plot and return best result
plt(tipBest.crds', 'mo', 10);
tipIdx  = tipBest.idx;
tipCrds = tipBest.crds;

end
