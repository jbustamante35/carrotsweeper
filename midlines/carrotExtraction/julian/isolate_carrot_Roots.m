function out = isolate_carrot_Roots(msk, vis, dOut, idx)
%% isolate_carrot_Roots: customized version of PhytoMorph's isolate_Roots
% This function is difficult to understand...
%
% Usage:
%   out = isolate_carrot_Roots(msk, vis, dOut, idx)
%
% Input:
%   msk: binary mask image
%   vis: boolean to visualize output
%   dOut: directory to output data
%   idx: array index of input image
%
% Output:
%   out: output structure containing midline and contour data
%

%% Constants and Data structure set-up
DEBUG = 0;

% Contstants for finding highest curvature
KSNIP         = 30;  % Original 50
SMOOTH_VALUE  = 100; % Original 15
SMOOTH_VALUE2 = 30; % Original 30

% Thresholds for filtering out small sizes of data
CURV_THRESH = 80;  % Original 80
SKEL_THRESH = 300; % Original 500
MASK_THRESH = 300; % Original 300

% Constants for tracking point at gradient
MAX_STEP = 50000;    % Original 50000
RHO      = 20;       % Original 20
RAD      = pi / 2;   % Original pi / 2
PDENSITY = [20 200]; % Original [20 200]
WSIGMA   = 0.3;      % Original 0.3

% Output data structures
curve    = struct('level', [], 'data', [], 'length', []);
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
    Io         = msk;
    
    %% Get skeleton structure and store contour
    % This is where everything goes wrong [ often creates empty data ]
    % Perhaps playing around with the SKEL_THRESH value will help?
    msk      = Io;
    gI       = double(imcomplement(bwdist(msk)));
    gI       = imfilter(gI,fspecial('gaussian', [31 31], 7));
    [g1, g2] = gradient(gI);
    thresh   = graythresh(msk);
    skel     = msk < thresh;
    skel     = bwareaopen(skel, SKEL_THRESH);
    
    % Setting edge with object to 0
    if sum(~msk(:, end)) > 0
        % Right side
        cIdx = size(msk, 2);
    else
        % Left sie
        cIdx = 1;
    end
    
    saveVec       = skel(:, cIdx);
    skel(:, cIdx) = 0;
    cB            = imclearborder(skel);
    cB(:, cIdx)   = saveVec;
    
    skel      = cB;
    skel      = imfill(skel, 'holes');
    
    %
    cB   = imclearborder(skel);
    skel = skel - cB;
    skel = imclose(skel, strel('disk', 21, 0));
    skel = imfill(skel, 'holes');
    
    C = contourc(double(skel), [1 1]);
    
    %% Get the curve structure
    str     = 1;
    c       = 1;
    ttlCrds = size(C, 2);
    
    while str < ttlCrds
        ed              = str + C(2, str);
        curve(c).level  = C(1, str);
        curve(c).data   = C(:, str+1:ed);
        curve(c).length = size(curve(c).data, 2);
        c               = c + 1;
        str             = ed + 1;
    end
    
    % Remove curves below threshold
    %     curve(curve.length < CURV_THRESH) = [];
    curve(cell2mat(arrayfun(@(x) x.length < CURV_THRESH, ...
        curve, 'UniformOutput', 0))) = [];
    
    %% Find the highest curvature
    numCurves = numel(curve);
    tipIdx    = cell(1, numCurves);
    for i = 1 : numCurves
        o                = cwtK(curve(i).data', {SMOOTH_VALUE});
        [~, tipIdx{i}]   = min(o.K);
        o                = cwtK(curve(i).data', {SMOOTH_VALUE2});
        [~, fine_tipIDX] = min((o.K(tipIdx{i} - KSNIP : tipIdx{i} + KSNIP)));
        tipIdx{i}        = tipIdx{i} + (fine_tipIDX - KSNIP - 1);
    end
    
    img = double(bwdist(~skel));
    g1  = -g1;
    g2  = -g2;
    for i = 1: numCurves
        
        tipCrds = curve(i).data(:, tipIdx{i});
        tipX    = tipCrds(1);
        tipY    = tipCrds(2);
        t1      = ba_interp2(g1, tipX, tipY);
        t2      = ba_interp2(g2, tipX, tipY);
        
        %
        N     = [t2 t1];
        N     = -N / norm(N);
        T     = [N(2) -N(1)];
        iDirc = [T ; N];
        
        %%
        %         quiver(curve(e).data(1,:),curve(e).data(2,:),d1X1,d1X2)
        %         quiver(curve(e).data(1,:),curve(e).data(2,:),-d1X2,d1X1,'r')
        
        midlines(i).data = trackFromPointAtGradient_carrot(img, tipCrds, iDirc, ...
            MAX_STEP, RHO, RAD, PDENSITY, WSIGMA);
        
    end
    
    %% Debug curvature finder with differing values of constants
    % A main point of failure is when the function to compute the point of
    % highest curvature identifies a coordinate on the side of the object,
    % rather than the desired tip of the structure. This debug function iterates
    % through a broad range of values for the following constants and plots
    % the highest curvature point. 
    %
    if DEBUG
        cla;clf;
        imshow(msk, []);
        hold on;
        rng               = 1 : 20 : 201;
        groundTruth       = [327 ; 234];
        [tipIdx, tipCrds] = runDebugMode(curve, rng, groundTruth, numCurves);
        
        %%
        tipX    = tipCrds(1);
        tipY    = tipCrds(2);
        img     = double(bwdist(~skel));
        g1      = -g1;
        g2      = -g2;        
        
        %
        t1    = ba_interp2(g1, tipX, tipY);
        t2    = ba_interp2(g2, tipX, tipY);
        N     = [t2 t1];
        N     = -N / norm(N);
        T     = [N(2) -N(1)];
        iDirc = [T ; N];
        
        %%
        midlines(i).data = trackFromPointAtGradient_carrot(img, tipCrds, iDirc, ...
            MAX_STEP, RHO, RAD, PDENSITY, WSIGMA);
    end
    
    %% Visualize outputs
    if vis
        hold off;
        image(cat(3,double(Io) / max(Io(:)), double(Io) / max(Io(:)), ...
            double(Io) / max(Io(:))));
        colormap('gray')
        axis off;
        hold on;
        
        for i = 1 : numCurves
            plot(curve(i).data(1,:),    curve(i).data(2,:),    'r');
            plot(midlines(i).data(1,:), midlines(i).data(2,:), 'g');
        end
        
        for i = 1 : numCurves
            plot(curve(i).data(1,tipIdx{i}), curve(i).data(2,tipIdx{i}), 'g*');
        end
        
        drawnow;
        
        if nargin == 4
            if ~isempty(dOut)
                mkdir(dOut);
                saveas(gca,[dOut idx '.tif']);
            end
        end
    end
    
    out.midlines = midlines;
    out.contours = curve;
    
catch e
    fprintf(2, 'Error isolating root\n%s\n', e.getReport);
    out.ME = e;
end

end

function [tipIdx, tipCrds] = runDebugMode(curve, rng, groundTruth, n)
%%
for k = rng
    col   = 'mo';
    err   = 1;
    for kk = rng
        for kkk = rng
            KSNIP         = k;   % Original 50
            SMOOTH_VALUE  = kk;  % Original 15
            SMOOTH_VALUE2 = kkk; % Original 30
            try
                numCurves = numel(curve);
                tipIdx    = cell(1, numCurves);
                o                = cwtK(curve(n).data', {SMOOTH_VALUE});
                [~, tipIdx{n}]   = min(o.K);
                o                = cwtK(curve(n).data', {SMOOTH_VALUE2});
                [~, fine_tipIDX] = min((o.K(tipIdx{n} - KSNIP : tipIdx{n} + KSNIP)));
                tipIdx{n}        = tipIdx{n} + (fine_tipIDX - KSNIP - 1);
                tipCrds          = curve(n).data(:, tipIdx{n});
            catch 
                tipCrds = groundTruth;
                err     = 2;
                col     = 'yo';
            end
            plt(tipCrds', col, 8);
            drawnow;
            fprintf(err, 'KSNIP %d | SMOOTH_VALUE %d | SMOOTH_VALUE2 %d\n', ...
                KSNIP, SMOOTH_VALUE, SMOOTH_VALUE2);
        end
    end
end
end