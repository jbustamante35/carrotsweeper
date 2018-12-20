function [DomainS, DomainG] = extendCarrotMidline(mline, dTrans, img, LP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub-sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate curvilinear-domain
WIDTH_NUMP = 200;
PCA_RHO = 15;
WIDTH = 200;
Domain = genCurvilinearDomain(mline, PCA_RHO, WIDTH, WIDTH_NUMP, img);

% extension on one side
dX = -diff(Domain,1,1);
SNIP = 50;
dX = mean(dX(1:SNIP,:,:),1);
dNOR = sum(dX.^2,3).^-.5;
dX = bsxfun(@times,dX,dNOR);

if nargin == 3
    LP = mean(sum(diff(mline,1,1).^2,2).^.5);
end

EXT = 20;
EXT = linspace(0,EXT,EXT/LP);
addedP = numel(EXT);
EXT = bsxfun(@times,EXT',dX);
EXT = bsxfun(@plus,EXT,Domain(1,:,:));
Domain = cat(1,flipdim(EXT,1),Domain);
% extension on one side
Domain = flipdim(Domain,1);
dX = -diff(Domain,1,1);
SNIP = 50;
dX = mean(dX(1:SNIP,:,:),1);
dNOR = sum(dX.^2,3).^-.5;
dX = bsxfun(@times,dX,dNOR);

if nargin == 3
    LP = mean(sum(diff(mline,1,1).^2,2).^.5);
end

EXT = 20;
EXT = linspace(0,EXT,EXT/LP);
addedP = numel(EXT);
EXT = bsxfun(@times,EXT',dX);
EXT = bsxfun(@plus,EXT,Domain(1,:,:));
Domain = cat(1,flipdim(EXT,1),Domain);
Domain = flipdim(Domain,1);

% reshape for sub-sampling
dsz = size(Domain);
DomainG = Domain;
DomainG(:,:,1) = DomainG(:,:,1) + dTrans(1);
DomainG(:,:,2) = DomainG(:,:,2) + dTrans(2);
DomainS = reshape(Domain,[dsz(1)*dsz(2) dsz(3)]);
DomainS = bsxfun(@plus,DomainS,dTrans);
end
