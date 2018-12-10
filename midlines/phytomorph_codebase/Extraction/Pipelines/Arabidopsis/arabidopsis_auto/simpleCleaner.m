function [D L] = simpleCleaner(D,L,absA,initA)
    ridx = find(any(D == 0,1) | any(abs(diff(D,1,1)) > absA*pi/180,1) | abs(D(1,:)) > initA*pi/180);
    D(:,ridx) = [];
    L(ridx) = [];
 %   D = bsxfun(@minus,D,D(1,:));
    %D = gradient(D')';
end