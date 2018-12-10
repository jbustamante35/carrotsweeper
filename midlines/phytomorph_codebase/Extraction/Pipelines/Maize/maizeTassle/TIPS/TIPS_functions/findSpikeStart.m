function [ firstBranch ] = findSpikeStart( width, tBin, spike, tol )
%FINDSPIKESTART Locates the lowest branchpoint on the tassel
%   Takes slices across within a window of specified size centered on the
%   spike in a binary image.  Sums of the slices are differentiated and the
%   point where the derivative exceeds 'tol' is taken as the first branch
%   point.
%       width: width of the windows centered on the spike
%       tBin: tassel binary image
%       spike: 2 x n matrix of pixels located on the spike
%       tol: tolerance for derivative of window sums.  Used to identify the
%            first branch.

%[~, x] = find(tBin);
bottom = min(spike(:,1));

if mod(width, 2) == 0
    width = width + 1;
end

halfwidth = (width - 1) / 2;

length = size(spike, 1);
slice = zeros(length, 1);
start = size(spike, 1);
j = 1;
for i=start:-1:1
    winStart = spike(i,2) - halfwidth;
    winEnd = spike(i,2) + halfwidth;

    slice(j,1) = sum(tBin(winStart:winEnd, spike(i,1)));
    j = j + 1;
end

sigma = 5;
range = -20:20;
K=inline('exp(-(x.^2)/2/sig^2)');
kern = K(sigma,range) / sum(K(sigma, range)); 

smooth = conv(slice, kern, 'same');

dsmooth = diff(smooth);
dsmooth = dsmooth((max(range)+1):size(dsmooth, 1));

firstBranch = find(dsmooth > tol, 1);
% Add two to compensate for the extra 1 in line 18 and non-overlap of
% kernel end and firstBranch start.
firstBranch = firstBranch + max(range) + bottom + 2;

end













