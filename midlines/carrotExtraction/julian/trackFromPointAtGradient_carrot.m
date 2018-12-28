function curve = trackFromPointAtGradient_carrot(skel, tipCrds, initDirection, MAX_STEP, RHO, RAD, PDENSITY, WSIGMA)
%% trackFromPointAtGradient_carrot: follow curve on path
% This is where I'd describe this function but I have no idea how it works
%
% Usage:
%   curve = trackFromPointAtGradient_carrot(I, P, initD, maxStep, RHO, RAD, pointDensity, wsigma)
%
% Input:
%   skel: processed skeleton image after thresholding
%   tip: x-/y-coordinates of object's tip
%   initDirection vector designating direction to begin tracking
%   MAX_STEP:
%   RHO:
%   RAD:
%   PDENSITY:
%   WSIGMA:
%
% Output:
%   curve:
%

%% Set-up for algorithm
x = tipCrds(1);
y = tipCrds(2);

%% Constants for Belief Function
% the width of the belief or the momentum is a function of percent of cutoff
MAX_WIDTH   = 60 * pi / 180;
K           = 1.5;
ALPHA       = 1;
SCALE       = 10;
STOP_THRESH = 300;

%% Function handles for Width, Radial, Angula Belief, and then Stop Function
width  = @(delta)           MAX_WIDTH * normpdf(delta, 0, WSIGMA) / normpdf(0, 0, WSIGMA);
radial = @(x)               wblpdf(x * SCALE^-1, ALPHA, K);
angle  = @(x, delta)        normpdf(x, 0, width(delta));
blief  = @(rad, rho, delta) angle(rad, delta) .* radial(rho);
stopFn = @(path)            size(path, 2) < MAX_STEP & path(1, end) > STOP_THRESH;

%% Create object containing curve data
T = goT();
T.setWfunction(blief);
T.setNhoodRho(RHO);
T.setNhoodRad(RAD);
T.setNhoodDensity(PDENSITY);
T.generateH();
T.setPosition([x;y]);
T.setImage(skel);
T.setDirection(initDirection);

% Follow curve until stop function
T.walkUntil(stopFn);
T.reparameterizeCurve();

curve = T.position;

end