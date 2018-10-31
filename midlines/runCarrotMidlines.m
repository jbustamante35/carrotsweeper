function M = runCarrotMidlines(directories)
%% runCarrotMidlines: run multiple directories through carrotMidlinesPipeline
% 
% 
% Usage:
% 
% 
% Input:
% 
% 
% Output:
% 
% 

D      = dir(directories);
D(1:2) = [];
D      = D(cat(1, D.isdir));
M      = arrayfun(@(x) carrotMidlinesPipeline([x.folder '/' x.name]), D, 'UniformOutput', 0);

end