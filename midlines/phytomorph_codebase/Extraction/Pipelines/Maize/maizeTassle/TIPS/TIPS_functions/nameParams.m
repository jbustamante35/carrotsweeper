function [ P ] = nameParams( params )
%NAMEPARAMS Assigns parameters to named variables
%   gThresh: threshold value below which to abort further analysis.  
P.('gThresh') = params(1);
P.('padSize') = params(2);
P.('smoothSigma') = params(3);
P.('smoothKernelDim') = params(4);
P.('skelTol') = params(5);
P.('skelMinBranch') = params(6);
P.('spikeWidth') = params(7);
P.('spikeTol') = params(8);

end

