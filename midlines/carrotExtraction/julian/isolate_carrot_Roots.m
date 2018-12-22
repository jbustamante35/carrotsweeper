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
% Contstants for finding highest curvature
KSNIP        = 50;
SMOOTH_VALUE = 15;
bint2        = @(g, cX, cY) ba_interp2(g, cX, cY);

% Thresholds for filtering out small sizes of data
CURV_THRESH = 100;
SKEL_THRESH = 500;
MASK_THRESH = 300;

% Constants for tracking point at gradient
MAX_STEP = 50000;
RHO      = 20;
RAD      = pi / 2;
PDENSITY = [20 200];
WSIGMA   = 0.3;

% Output data structures
curve    = struct('level', [], 'data', [], 'length', []);
midlines = struct('data', []);
out      = struct('midlines', [], 'contours', []);

try
    %% Pre-process mask by facing object left-right
    if size(msk, 3) > 1
        msk = rgb2gray(msk);
    end
    
    msk(:,end) = [];
    msk        = handleFLIP(msk, []);
    msk        = padarray(msk, [0 MASK_THRESH], 'pre', 'replicate');
    Io         = msk;
    
    %% Get skeleton structure and store contour
    % This is where everything goes wrong [ often creates empty data ]
    % Perhaps playing around with the SKEL_THRESH value will help?
    gI       = double(imcomplement(bwdist(msk)));
    gI       = imfilter(gI,fspecial('gaussian', [31 31], 7));
    [g1, g2] = gradient(gI);
    thresh   = graythresh(msk);
    skel     = msk < thresh;
    skel     = bwareaopen(skel, SKEL_THRESH);
    
    %
    saveVec   = skel(:, 1);
    skel(:,1) = 0;
    cB        = imclearborder(skel);
    cB(:,1)   = saveVec;
    skel      = cB;
    skel      = imfill(skel,'holes');
    
    %
    cB   = imclearborder(skel);
    skel = skel - cB;
    skel = imclose(skel,strel('disk', 21, 0));
    skel = imfill(skel, 'holes');
    
    C = contourc(double(skel), [1 1]);
    
    %% Get the curve structure
    str       = 1;
    c         = 1;
    ttlCurves = size(C, 2);
    
    while str < ttlCurves
        ed              = str + C(2, str);
        curve(c).level  = C(1, str);
        curve(c).data   = C(:, str+1:ed);
        curve(c).length = size(curve(c).data, 2);
        c               = c + 1;
        str             = ed + 1;
    end
    
    % Remove curves below threshold
    curve(curve.length < CURV_THRESH) = [];
    
    %% Find the highest curvature
    numCurves = numel(curve);
    tipIDX    = cell(1, numCurves);
    for e = 1 : numCurves
        o                = cwtK(curve(e).data', {SMOOTH_VALUE});
        [~, tipIDX{e}]   = min(o.K);
        o                = cwtK(curve(e).data', {30});
        [~, fine_tipIDX] = min((o.K(tipIDX{e} - KSNIP : tipIDX{e} + KSNIP)));
        tipIDX{e}        = tipIDX{e} + (fine_tipIDX - KSNIP - 1);
    end
    
    skelImg = double(bwdist(~skel));
    g1      = -g1;
    g2      = -g2;
    for e = 1: numCurves
        
        %
        crds  = curve(e).data(:, tipIDX{e});
        crdsX = curve(e).data(1, tipIDX{e});
        crdsY = curve(e).data(2, tipIDX{e});
        t1    = bint2(g1, crdsX, crdsY);
        t2    = bint2(g2, crdsX, crdsY);
        
        %
        N     = [t2 t1];
        N     = -N / norm(N);
        T     = [N(2) -N(1)];
        initD = [T ; N];
        
        %
        %quiver(curve(e).data(1,:),curve(e).data(2,:),d1X1,d1X2)
        %quiver(curve(e).data(1,:),curve(e).data(2,:),-d1X2,d1X1,'r')
        
        midlines(e).data = trackFromPointAtGradient_carrot(skelImg, crds, initD, ...
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
        
        for e = 1 : numCurves
            plot(curve(e).data(1,:),    curve(e).data(2,:),    'r');
            plot(midlines(e).data(1,:), midlines(e).data(2,:), 'g');
        end
        
        for e = 1 : numCurves
            plot(curve(e).data(1,tipIDX{e}), curve(e).data(2,tipIDX{e}), 'g*');
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
    
catch err
    fprintf('Error isolating root\n%s\n', err.getReport);
    out.ME = err;
end

end