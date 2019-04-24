function [domainS, domainG] = extendCarrotMidline(mline, domainTranslation, msk, np, lp)
%% extendCarrotMidline:
% I think this function extends from the midline to the contour boundary?
%
% Usage:
%
%
% Input:
%   mline: coordinates of the path along a midline
%   domainTranslation:
%   msk: binary mask image for corresponding midline
%   wid: width to extend interpolation of straightened mask
%   lp:
%
% Output:
%   domainS:
%   domainG:
%

%% Generate curvilinear-domain
LPARGS     = 4; % number of arguments if lp is used
% WIDTH_NUMP = 200; % Dynamically set to original mask width
PCA_RHO    = 15;
% WIDTH      = 200; % Dynamically set to original mask width
wid = round(np / 2); % 
% domain     = genCurvilinearDomain(mline, PCA_RHO, WIDTH, WIDTH_NUMP, msk, 0);
domain     = genCurvilinearDomain(mline, PCA_RHO, wid, np, msk, 0);

%% extension on one side
dX   = -diff(domain, 1, 1);
SNIP = 50;
dX   = mean(dX(1:SNIP,:,:), 1);
dNOR = sum(dX.^2, 3).^-0.5;
dX   = bsxfun(@times, dX, dNOR);

if nargin == LPARGS
    lp = mean(sum(diff(mline, 1, 1).^2, 2).^0.5);
end

%%
EXT    = 20;
EXT    = linspace(0, EXT, EXT / lp);
EXT    = bsxfun(@times, EXT', dX);
EXT    = bsxfun(@plus, EXT, domain(1,:,:));
domain = cat(1, flip(EXT, 1), domain);

%% extension on one side
domain = flip(domain, 1);
dX     = -diff(domain, 1, 1);
SNIP   = 50;
dX     = mean(dX(1:SNIP,:,:), 1);
dNOR   = sum(dX.^2,3).^-.5;
dX     = bsxfun(@times, dX, dNOR);

if nargin == LPARGS
    lp = mean(sum(diff(mline, 1, 1).^2, 2).^0.5);
end

%%
EXT    = 20;
EXT    = linspace(0, EXT, EXT / lp);
EXT    = bsxfun(@times, EXT', dX);
EXT    = bsxfun(@plus, EXT, domain(1,:,:));
domain = cat(1, flip(EXT,1), domain);
domain = flip(domain,1);

%% Reshape for sub-sampling
dsz            = size(domain);
domainG        = domain;
domainG(:,:,1) = domainG(:,:,1) + domainTranslation(1);
domainG(:,:,2) = domainG(:,:,2) + domainTranslation(2);
domainS        = reshape(domain,[dsz(1) * dsz(2) dsz(3)]);
domainS        = bsxfun(@plus, domainS, domainTranslation);

end