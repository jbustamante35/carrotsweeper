clear all

Xi = linspace(0,2,100);
syms x xo k ix fx
df(x,xo,k,ix,fx) = (fx-ix).*(1+exp(-k*(x-xo))).^-1 + ix;

%%
pulseVec = sym('d', [1 4]);
notVec = [-1 1 1.4];
rateVec = [50000 100 100];
fout = createTrain(df,pulseVec,notVec,rateVec);
%%

Xi = linspace(0,2,100);
pV = [0 8 20 30];
pV = [0 diff(pV)];
Yi = fout(Xi,pV(1),pV(2),pV(3),pV(4));
close all
for e = 1:size(Yi)
    hold on
    plot(Xi,Yi{e})
end
F = sum(fout);
plot(Xi,F(Xi,pV(1),pV(2),pV(3),pV(4)));
sfout = sum(fout);
%%
ifout = int(sum(fout),x);
iY = ifout(Xi,pV(1),pV(2),pV(3),pV(4)) - ifout(0,pV(1),pV(2),pV(3),pV(4));
close all
plot(Xi,iY)
%%
syms solveRate;
max = 40;
xRate = [0 8 20 50];
yRate = [0 16 20 50];
xRate = [0 diff(xRate)];
yRate = [0 diff(yRate)];
solveAt = 1.4;
X = ifout(solveAt,xRate(1),xRate(2),xRate(3)+solveRate,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+solveRate,xRate(4));
Y = ifout(solveAt,yRate(1),yRate(2),yRate(3)-solveRate,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-solveRate,yRate(4));
preX = ifout(Xi,xRate(1),xRate(2),xRate(3),xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3),xRate(4));
preY = ifout(Xi,yRate(1),yRate(2),yRate(3),yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3),yRate(4));
preRateX = sfout(Xi,xRate(1),xRate(2),xRate(3),xRate(4));
preRateY = sfout(Xi,yRate(1),yRate(2),yRate(3),yRate(4));
%%
for e = 1:10
    
    solveAt = 1.4;
    X = ifout(solveAt,xRate(1),xRate(2),xRate(3)+solveRate,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+solveRate,xRate(4));
    Y = ifout(solveAt,yRate(1),yRate(2),yRate(3)-solveRate,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-solveRate,yRate(4));

    r = solve(X==Y,solveRate);
    nextPointX = vpa(ifout(notVec(3),xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4)));
    nextPointY = vpa(ifout(notVec(3),yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4)));
    nextPointX
    nextPointY

    solveAt = 2;
    X = ifout(solveAt,xRate(1),xRate(2),xRate(3)+r,xRate(4)+solveRate) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4)+solveRate);
    Y = ifout(solveAt,yRate(1),yRate(2),yRate(3)-r,yRate(4)-solveRate) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4)-solveRate);
    r2 = solve(X==Y,solveRate)
    xRate(4) = xRate(4) + r2;
    yRate(4) = yRate(4) - r2;
end
%%
Xs = ifout(Xi,xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4));
Ys = ifout(Xi,yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4));
postX = ifout(Xi,xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4));
postY = ifout(Xi,yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4));
close all
plot(Xi,preX,'c')
hold on
plot(Xi,preY,'k')
plot(Xi,postX,'b--')
plot(Xi,postY,'g--')


postRateX = sfout(Xi,xRate(1),xRate(2),xRate(3),xRate(4));
postRateY = sfout(Xi,yRate(1),yRate(2),yRate(3),yRate(4));
figure;
plot(Xi,preRateX,'k')
hold on
plot(Xi,postRateX,'r--')
plot(Xi,preRateY,'c')
plot(Xi,postRateY,'b--')
%%

Baseline = 10;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline+1;

clear X
[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,2),X(:,:,1));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;
rhoX = ifout(G(:,:,2),xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4));
rhoY = ifout(G(:,:,2),yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4));

rhoX = double(vpa(rhoX));
rhoY = double(vpa(rhoY));

rhoLine = linspace(1.4,2,100);
rhoConst = 1.4*ones(size(rhoLine));


wowRho = (rhoLine.^2 + rhoConst.^2).^.5;
wowRho = 1.8*ones(size(rhoLine));


whatRhoX = ifout(wowRho,xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4));
whatRhoY = ifout(wowRho,yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4));

whereTH = atan2(rhoLine,rhoConst);
whereTH = linspace(0,pi/2,100);
whyRhoX = double(vpa(whatRhoX.*cos(whereTH)));
whyRhoY = double(vpa(whatRhoY.*sin(whereTH)));


rhoXde = ifout(linspace(1.4,2,10),xRate(1),xRate(2),xRate(3)+r,xRate(4)) - ifout(0,xRate(1),xRate(2),xRate(3)+r,xRate(4));
rhoYde = ifout(linspace(1.4,2,10),yRate(1),yRate(2),yRate(3)-r,yRate(4)) - ifout(0,yRate(1),yRate(2),yRate(3)-r,yRate(4));

rhoXde = double(vpa(rhoXde));
rhoYde = double(vpa(rhoYde));

clear nG
nG(:,:,1) = rhoX.*cos(G(:,:,1));
nG(:,:,2) = rhoY.*sin(G(:,:,1));


V = ba_interp2(double(nG(:,:,1)),whyRhoY,whyRhoX)
close all
CL = {'r'};
skip = 1;
for e1 = 1:skip(1):size(nG,2)
    plot(nG(:,e1,1),nG(:,e1,2),CL{1})
    hold on
end
plot(nG(:,end,1),nG(:,end,2),CL{1})

skip = [1];
for e1 = 1:skip(1):size(nG,1)
    plot(nG(e1,:,1),nG(e1,:,2),CL{1})
    hold on
end
plot(nG(end,:,1),nG(end,:,2),CL{1})

plot(whyRhoX,whyRhoY,'b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
syms x xo0 k0 ix0 fx0
f0 = df(x,xo0,k0,ix0,fx0);
syms x xo1 k1 ix1 fx1
f1 = df(x,xo1,k1,ix1,fx1);
d = f0 + f1;
%%


[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,2),X(:,:,1));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;

%%
%df = taylor(df,x,'Order',50);
f = int(df,x);
ddf = diff(df,x);


syms x transitiondF1 k1 initSlope1 finalSlope1
f1(x,transitiondF1,k1,initSlope1,finalSlope1) = df(x,transitiondF1,k1,initSlope1,finalSlope1);
syms x transitiondF2 k2 initSlope2 finalSlope2
f2(x,transitiondF2,k2,initSlope2,finalSlope2) = df(x,transitiondF2,k2,initSlope2,finalSlope2);
syms blendPoint blendRate v1 v2
blenddF(x,blendPoint,blendRate,v1,v2) = v1*df(x,blendPoint,blendRate,1,0) + v2*df(x,blendPoint,blendRate,0,1);



in_f1 = argnames(f1);
in_f2 = argnames(f2);
in_blend = argnames(blenddF);
in_dF = [in_blend(:);in_f1(:);in_f2(:)];

dF = myCompose(blenddF,f1,v1,[],[]);
dF = myCompose(dF,f2,v2,[],[]);



xValue = Xi;
blendPointValue = 1.4;
blendRateValue = 50;

transitiondF1Value = 1;
k1Value = 50;
initSlope1Value = 8;
finalSlope1Value = 22;

transitiondF2Value = 1;
k2Value = 50;
initSlope2Value = 0;
finalSlope2Value = 16;


dr = dF(xValue,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value);
R = int(dF,x);

%%
R = int(dF,x,0,x);
nm = argnames(R);
R([x;nm(:)]) = R;
%%
%{
%%
R = tmp(xValue,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value) - ...
  tmp(0,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value);
%%
R = tmp;
%}
%%

fargs = argnames(dF);
dF = matlabFunction(dF,'Vars',fargs);

close all
plot(Xi,double(dr),'r');
hold on
plot(Xi,double(dr),'g--');
F = @(d,xValue,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value)...
  cumsum(dF(xValue,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value),d);
Y = F(2,xValue,blendPointValue,blendRateValue,...
  transitiondF1Value,k1Value,initSlope1Value,finalSlope1Value,...
  transitiondF2Value,k2Value,initSlope2Value,finalSlope2Value);
figure;
plot(Xi,Y,'b--');
%%

%%
smallMax = 2;
Xi = linspace(0,smallMax,100);

Xlim = 50;
Ylim = Xlim;

XblendPointValue = 1.4;
XblendRateValue = 50;

XtransitiondF1Value = 1;
Xk1Value = 50;
XinitSlope1Value = 12;
XfinalSlope1Value = 22;

XtransitiondF2Value = 1;
Xk2Value = 50;
XinitSlope2Value = 0;
XfinalSlope2Value = 16;

YblendPointValue = 1.4;
YblendPointValue = XblendPointValue;
YblendRateValue = 50;

YtransitiondF1Value = 1;
Yk1Value = 50;
YinitSlope1Value = 8;
YfinalSlope1Value = 22;

YtransitiondF2Value = 1;
Yk2Value = 50;
YinitSlope2Value = 0;
YfinalSlope2Value = 16;

SolveAt = XblendPointValue;
matchedFinalSlope = (smallMax - XblendPointValue)*Xlim;
YfinalSlope2Value = matchedFinalSlope;
XfinalSlope2Value = matchedFinalSlope;

probeX = R(Xi,XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value);
probeX = double(vpa(eval(probeX)));


probeY = R(Xi,XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value);
probeY = double(vpa(eval(probeY)));

DprobeX = dF(Xi,XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value);
DprobeX = double(vpa((DprobeX)));

DprobeY = dF(Xi,YblendPointValue,YblendRateValue,...
    YtransitiondF1Value,Yk1Value,YinitSlope1Value,YfinalSlope1Value,...
    YtransitiondF2Value,Yk2Value,YinitSlope2Value,YfinalSlope2Value);
DprobeY = double(vpa((DprobeY)));

close all
plot(Xi,probeX,'r');
hold on
plot(Xi,probeY,'g--');
figure;

plot(Xi,DprobeX,'r');
hold on
plot(Xi,DprobeY,'g--');

syms XexpandSlope;
EXP = R(SolveAt,XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value) == ...
R(SolveAt,YblendPointValue,YblendRateValue,...
    YtransitiondF1Value,Yk1Value,YinitSlope1Value,XexpandSlope,...
    YtransitiondF2Value,Xk2Value,YinitSlope2Value,YfinalSlope2Value);

YexpandSlope = vpa(simplify(solve(EXP,XexpandSlope)));
double(vpa(eval(EXP)))
%%

%%


Major = F(1,Xi,XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value);
close all
plot(Major)
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Baseline = 40;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline+1;



[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
P(:,:,2) = atan2(X(:,:,2),X(:,:,1));
P(:,:,1) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;

M = F(1,P(:,:,1),XblendPointValue,XblendRateValue,...
    XtransitiondF1Value,Xk1Value,XinitSlope1Value,XfinalSlope1Value,...
    XtransitiondF2Value,Xk2Value,XinitSlope2Value,XfinalSlope2Value);
%%

m = F(P(:,:,1),YblendPointValue,YblendRateValue,...
    YtransitiondF1Value,Yk1Value,YinitSlope1Value,YfinalSlope1Value,...
    YtransitiondF2Value,Yk2Value,YinitSlope2Value,YfinalSlope2Value);
m = simplify(m);
%%
clear all
syms theta R
x(theta,R) = piecewise(theta >= -pi/4 & theta <= pi/4,R,...
                        theta > pi/4 & theta < 3*pi/4,R-R*sin(theta),...
                        theta > -3*pi/4 & theta < -pi/4,R-R*sin(theta),...
                        theta <= pi & theta >= 3*pi/4 | theta >= pi & theta <= -3*pi/4,-R);
y(theta,R) = piecewise(theta >= -pi/4 & theta <= pi/4,R-R*cos(theta),...
                        theta > pi/4 & theta < 3*pi/4,R,...
                        theta > -3*pi/4 & theta < -pi/4,-R,...
                        theta <= pi & theta >= 3*pi/4 | theta >= pi & theta <= -3*pi/4,R-R*cos(theta));
S(theta,R) = [x(theta,R);y(theta,R)];


%%
close all
TH = linspace(-pi/4,pi/4,100);
for e = 1:numel(TH)
    F(:,e) = double(vpa(S(TH(e),4)));
end
plot(F(1,:),F(2,:),'.')
hold on

TH = linspace(pi/4,3*pi/4,100);
for e = 1:numel(TH)
    F(:,e) = double(vpa(S(TH(e),4)));
end
plot(F(1,:),F(2,:),'.')













