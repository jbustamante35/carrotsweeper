function [ nb, C, lowBranch, r ] = BNcircle( tBin, spike, branchpoints, startpoint, firstBranch )
%BNCIRCLE counts branch number by creating concentric circles around a
%given point and counting the number of times they cross the binary tassel
%image.
%   tBin: binary tassel image
%   spike: coordinate list of pixels along the skeleton that comprise the
%      central spike
%   branchpoints: coordinates of all skeleton branch points
%   startpoint: if not using the lowest branchpoint as the origin, you can
%      specify whether to use the base or center of mass as the origin
%      instead ('base', 'center', respectively).

if ~isempty(branchpoints)
    %D = pdist2(branchpoints, spike);
    %D = min(D, [], 2);
    %branches = branchpoints(D <= 2, :);
    
    %[~, lowBranch] = min(branches(:,1));
    %lowBranch = branches(lowBranch,:); 
    
    % Try starting from point used for TLadj
    lowBranch = firstBranch;
    lowBranch = spike(lowBranch, :);
    
elseif isempty(branchpoints) && isempty(startpoint)
    % If the tassel has no branchpoints but a startpoint wasn't specified,
    % default to the base
    [~, lowBranch] = min(spike(:,1));
    lowBranch = spike(lowBranch,:);
elseif strcmp(startpoint, 'base')
    [~, lowBranch] = min(spike(:,1));
    lowBranch = spike(lowBranch,:);
elseif strcmp(startpoint,'center')
    lowBranch = regionprops(tBin, 'Centroid');
    lowBranch = lowBranch.Centroid;   
else
    msgID = 'BNcircle:badInputs';
    msg = 'unclear input options for BNcircle';
    exception = MException(msgID, msg);
    throw(exception);
end

edgeDist = abs([lowBranch(1) - [0 size(tBin, 2)] lowBranch(2) - [0 size(tBin, 1)]]);
edgeDist = max(edgeDist);

nb = [];

if numel(lowBranch) >= 1
    for r=100:50:edgeDist
    
        T = tBin;
    
        [x, y] = meshgrid(1:size(tBin,2), 1:size(tBin, 1));

        C = sqrt( (x-lowBranch(1)).^2 + (y - lowBranch(2)).^2 ) - r;  

        C = C <= 1 & C >= 0;

        T(C ~= 1) = 0;

        n = size(regionprops(T), 1);
        
        % Terminate loop if zero intersections with tassel
        if n == 0 && ~strcmp(startpoint, 'center'); break; end
    
        nb = [nb n];
    
    end
    
    [~, r] = max(nb);
    r = 100 + (r-1)*50;
    C = sqrt( (x-lowBranch(1)).^2 + (y - lowBranch(2)).^2 ) - r;  
    C = C <= 1 & C >= 0;
    [Cy, Cx] = find(C);
    C = [Cx Cy];

else
   nb = NaN;
   C = NaN;
   lowBranch = NaN;

end

