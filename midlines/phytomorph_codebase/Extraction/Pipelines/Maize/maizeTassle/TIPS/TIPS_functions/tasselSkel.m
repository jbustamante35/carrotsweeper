function [ endpoints, branchpoints, splines, lengths, maxPath, skelLength, base] = tasselSkel( smooth, tol, minBranchLength )
%TASSELSKEL skeletonizes the tassel image and fits splines from the base to
%all endpoints
%   smooth: smoothed binary tassel image
%   tol: tolerance - parameter for fitting splines
%   minBranchLength: minimum allowable length for branches.  Anything
%       shorter will be ignored.

skel = bwmorph(smooth, 'thin', Inf);
skelLength = sum(sum(skel));

adj = skel2adj(skel);

[y, x] = find(skel);

[endY, endX] = find(bwmorph(skel, 'endpoints'));
[~, endI] = ismember([endY, endX], [y, x], 'rows');
endI = endI(endI ~= 0);
endX = endX(endI ~= 0);
endY = endY(endI ~= 0);

endpoints = [endX, endY];
[bY, bX] = find(bwmorph(skel, 'branchpoints'));
branchpoints = [bX, bY];

% Find the base using an index that combines the distance from the
% bottom of the image and the distance from the center of the image
% This could use some fixing - won't always work

d2center = abs(y - (size(smooth, 1)/2));
baseIndx = x + d2center;

base = find(baseIndx == min(baseIndx), 1);



% For viewing endpoints, splines, etc:
%imshow(smooth); hold on;
%scatter(endX, endY, 100, 'b', '+')
%scatter(x(base,1), y(base,1), 100, 'r', '*')

% Distance from all skel points to closest branchpoint
D = bwdistgeodesic(skel, bX, bY); 

i = 1; % iterator for splines
keepE = false(size(endY, 1), 1); % indices of endpoints to keep (will exclude "short" branches)
maxPath = [];
for e = 1:size(endY, 1)
    
    if endI(e) == base
        continue % Skips making a spline from base to itself
    end
    
    [path, cost] = dijkstra(adj, base, endI(e));
    
    pathY = y(path);
    pathX = x(path);
    pathI = sub2ind(size(skel), y(path), x(path));
    
    if(size(pathY, 1) > size(maxPath, 1)); maxPath = [pathX, pathY]; end
    
    
    if size(pathI) < 50
        continue % skip if small distance to base - likely not a branch
    
    elseif D(endY(e), endX(e)) < minBranchLength 
        continue % also skip any branches with small distance to closest
    
    else
        % Parametrize the path and fit a spline
        toFit = [pathX, pathY];
        spline = csaps(1:size(toFit,1), toFit', tol);
        %fnplt(spline, 3, 'r');
    end
    
    if(i == 1)
        splines = spline;
        lengths = size(path, 1);
    else
        splines(i) = spline;
        lengths(i) = size(path, 1);
    end
    
    keepE(e) = 1;
    i = i + 1;
end

% If all endpoints are <75 pix from nearest branchpoint, no splines will be
% formed above.  This chunk fits a spline from the base to the endpoint
% farthest from it on the skeleton.
if ~exist('splines', 'var')
    D = bwdistgeodesic(skel, x(base,1), y(base,1));
    [farY farX] = find(D == max(max(D)));
    e = ismember([endX endY], [farX farY], 'rows');
    [path, cost] = dijkstra(adj, base, endI(e));
    pathY = y(path);
    pathX = x(path);
    pathI = sub2ind(size(skel), y(path), x(path));
    toFit = [pathX, pathY];
    spline = csaps(1:size(toFit,1), toFit', tol);
    splines = spline;
    lengths = size(path, 1);
end


% Return only endpoints from 
endpoints = endpoints(keepE, :);
% Return coordinates of base
base = [x(base, 1), y(base, 1)];

end

