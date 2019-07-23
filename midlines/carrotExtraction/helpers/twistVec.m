function tV = twistVec(vec)
%% twitstVec: twists 2D vector "in the plane" counter-clockwise
%
% Usage:
%   tV = twistVec(vec)
%
% Input:
%   vec: vector
%
% Output:
%   tV: twisted vector
%

%%
if size(vec,1) == 1
    tV = [-vec(2) , vec(1)];
    
elseif size(vec,2) == 1
    tV = [-vec(2) ; vec(1)];
end

%%
if numel(vec) == 3
    tV(end+1) = 0;
end
end
