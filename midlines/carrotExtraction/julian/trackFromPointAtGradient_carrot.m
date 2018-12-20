function curve = trackFromPointAtGradient_carrot(skel, crds, initD, maxStep, RHO, RAD, pointDensity, wsigma)
%% trackFromPointAtGradient_carrot: follow curve on path
% This is where I'd describe this function but I have no idea how it works
%
% Usage:
%   curve = trackFromPointAtGradient_carrot(I, P, initD, maxStep, RHO, RAD, pointDensity, wsigma)
%
% Input:
%   skel: processed skeleton image after thresholding
%   crds: x-/y-coordinates of curve to compute
%   initD: 
%   maxStep: 
%   RHO: 
%   RAD: 
%   pointDensity: 
%   wsigma: 
%
% Output:
%   curve:
%

%% Set-up for algorithm
x = crds(1);
y = crds(2);

%% Constants for Belief Function
% the width of the belief or the momentum is a function of percent of cutoff
maxWidth    = 60 * pi / 180;
k           = 1.5;
alpha       = 1;
scale       = 10;
STOP_THRESH = 300;

%% Function handles for Width, Radial, Angula Belief, and then Stop Function
width  = @(delta)           maxWidth * normpdf(delta, 0, wsigma) / normpdf(0, 0, wsigma);
radial = @(x)               wblpdf(x * scale^-1, alpha, k);
angle  = @(x, delta)        normpdf(x, 0, width(delta));
blief  = @(rad, rho, delta) angle(rad, delta) .* radial(rho);
stp    = @(path)            size(path, 2) < maxStep & path(1, end) > STOP_THRESH;

%% Create object containing curve data
T = goT();
T.setWfunction(blief);
T.setNhoodRho(RHO);
T.setNhoodRad(RAD);
T.setNhoodDensity(pointDensity);
T.generateH();
T.setPosition([x;y]);
T.setImage(skel);
T.setDirection(initD);

% Follow curve until stop function
T.walkUntil(stp);
T.reparameterizeCurve();

curve = T.position;

end