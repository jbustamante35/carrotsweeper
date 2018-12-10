function [f] = logistic(x,L,k,xo)
    f = L.*(1+exp(-k*(x-xo))).^-1;
end

%{

x = linspace(0,2,100);
f = logistic(x,2,30,1) + 8;
plot(x,f)


k = 10;

Baseline = 10;

STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline;


M = 14;
m = 8;



mMdelta = M - m;

D = @(th)mMdelta/2*cos(th*2)+mMdelta/2 - (Baseline - m);

th = linspace(-pi,pi,100);
plot(th,D(th));

mag1 = @(th)mMdelta/2*cos(th*2) + m + mMdelta/2;
mag2 = @(rho,d)d - logistic(rho,d,k,1) + Baseline;

rho = linspace(0,2,100);
plot(rho,mag2(rho,-2));

Q = linspace(0,2,100);
plot(Q,mag2(Q,2).*Q);

delta = pi/2;
MAX = pi;
y = @(x)mMdelta/2*abs(abs(x)-delta)/MAX+mMdelta/2 - (Baseline - m);
y = @(x)mMdelta/2*2*(abs(abs(x)-delta)/MAX*2 - .5)+mMdelta/2 - (Baseline - m);
D = y;
plot(linspace(-pi,pi,100),y(linspace(-pi,pi,100)));
hold on
plot(linspace(-pi,pi,100),D(linspace(-pi,pi,100)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
Baseline = 10;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline+1;


M = 14;
m = 8;


[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,2),X(:,:,1));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;

%%%


x = linspace(0,1,100)
Rx(2)-Rx(1)
v = i_logistic(x,20,100,1) + x; 
plot(x,v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all




Xi = linspace(0,2,100);
THi = linspace(0,2,100);

Baseline = 10;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline+1;

[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,1),X(:,:,2));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;

syms f dg x xi v w y init final rho theta f_rhoRate k i_rhoRate alpha rho_o i_minor_rhoRate i_major_rhoRate f_minor_rhoRate f_major_rhoRate

syms x xo0 xo1 kf0 kf1 kb0 kb1 if0 if1 ff0 ff1 fx xo ix xi0 xi1 

dg(rho,rho_o,k,i_rhoRate,f_rhoRate) = (f_rhoRate-i_rhoRate).*(1+exp(-k*((rho)-rho_o))).^-1 + i_rhoRate;

clear all
syms deltaSlope;
finalF = 30;
kMaster = 50;

k1 = kMaster;
initSlope1 = 8;
finalSlope1 = finalF - initSlope1;
transitionF1 = 1;

k2 = kMaster;
initSlope2 = 0;
finalSlope2 = finalSlope1;
transitionF2 = 1;


clear all
syms x xo k ix fx
df(x,xo,k,ix,fx) = (fx-ix).*(1+exp(-k*(x-xo))).^-1 + ix;

syms x transitionF1 k1 initSlope1 finalSlope1
f1(x,transitionF1,k1,initSlope1,finalSlope1) = df(x,transitionF1,k1,initSlope1,finalSlope1);
syms x transitionF2 k2 initSlope2 finalSlope2
f2(x,transitionF2,k2,initSlope2,finalSlope2) = df(x,transitionF2,k2,initSlope2,finalSlope2);
syms blendPoint blendRate v1 v2
blendF(x,blendPoint,blendRate,v1,v2) = v1*df(x,blendPoint,blendRate,1,0) + v2*df(x,blendPoint,blendRate,0,1);

in_f1 = argnames(f1);
in_f2 = argnames(f2);
in_blend = argnames(blendF);
in_F = [in_blend(:);in_f1(:);in_f2(:)];

syms l
F = myCompose(blendF,f1,v1,[],[]);
F = myCompose(F,f2,v2,[],[]);

xValue = Xi;
blendPointValue = 1.4;
blendRateValue = 50;

transitionF1Value = 1;
k1Value = 50;
initSlope1Value = 8;
finalSlope1Value = 22;

transitionF2Value = 1;
k2Value = 50;
initSlope2Value = 0;
finalSlope2 = 26;

F(xValue,blendPointValue,blendRateValue,...
    transitionF1Value,k1Value,initSlope1Value,finalSlope1Value,transitionF2Value


w = dg(Xi,1,100,.2,1);
v = dg(Xi,1,100,.2,1.2);


g = int(dg,rho);


g(rho,rho_o,k,i_rhoRate,f_rhoRate) = g - g(0,rho_o,k,i_rhoRate,f_rhoRate);
blend(v,w,x,xi,k) = v.*dg(x,xi,k,0,1) + w.*dg(x,xi,k,1,0);


plot(Xi,g(Xi,1,5,5,10));

k=30;
T = blend(v,w,Xi,1.4,30);


xPos(x,y,rho,theta,rho_o,k,i_minor_rhoRate,i_major_rhoRate,f_minor_rhoRate,f_major_rhoRate) = cos(theta).*g(rho,rho_o,k,i_minor_rhoRate,f_minor_rhoRate);
yPos(x,y,rho,theta,rho_o,k,i_minor_rhoRate,i_major_rhoRate,f_minor_rhoRate,f_major_rhoRate) = sin(theta).*g(rho,rho_o,k,i_major_rhoRate,f_major_rhoRate);


rot = [[cos(alpha);sin(alpha)],[-sin(alpha);cos(alpha)]];
pos(x,y,rho,theta,alpha,rho_o,k,i_minor_rhoRate,i_major_rhoRate,f_minor_rhoRate,f_major_rhoRate) = [xPos;yPos];%rot*[xPos;yPos];

test = g(Xi,1,20,10,20);
al = g(G(:,:,2),20,30,0,0);
%TY = xPos(X(:,:,2),X(:,:,1),G(:,:,2),G(:,:,1),1,k,m,M,F-m,F-M);

rhoFunc = matlabFunction(pos);

M = 14;
m = 8;
F = 60;
k = 200;
alpha = 0;


XX = rhoFunc(X(:,:,2),X(:,:,1),G(:,:,2),G(:,:,1),alpha,1,k,m,M,F-m,F-M);
XX = rhoFunc(X(:,:,2),X(:,:,1),G(:,:,2),G(:,:,1),alpha,1,k,m,M,F-m,F-M);

nG(:,:,1) = XX(1:41,:);
nG(:,:,2) = XX(42:end,:);


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



rhoFunc


yR = rhoFunc(G(:,:,2),1,k,m,MV-m);

w = solve(g(x,1,k,M,MV-M)==M,x);
xv = linspace(0,2,100);
yv = g(xv,1,10,10,20);
plot(xv,yv);

Baseline = 30;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline+1;

k = 55;
M = 14;
m = 8;
MV = 30;

rhoFunc = matlabFunction(g);
[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,2),X(:,:,1));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;


yR = rhoFunc(G(:,:,2),1,k,m,MV-m);
xR = rhoFunc(G(:,:,2),1,k,M,MV-M);


clear nG
nG(:,:,1) = xR.*cos(G(:,:,1));
nG(:,:,2) = yR.*sin(G(:,:,1));



%%%%%%%%%%%%%%%%%%%%%%%%
clear all
Baseline = 30;
STR = -2;
STP = 2;
STEP = (STP - STR)*Baseline;

syms f x L k xo alpha xLim
f(x) = .5*(x+x.^-1);
[X(:,:,1),X(:,:,2)] = ndgrid(linspace(STR,STP,STEP),linspace(STR,STP,STEP));
G(:,:,1) = atan2(X(:,:,2),X(:,:,1));
G(:,:,2) = (X(:,:,1).^2 + X(:,:,2).^2).^.5;
clear G
[G(:,:,2),G(:,:,1)] = ndgrid(linspace(1.2,2.3,50),linspace(-pi,pi,STEP));



Q = G(:,:,2).*exp(i*G(:,:,1));



TH = linspace(-pi,pi,100);
R = 1.2;
Q = R*exp(i*TH);
Q = f(Q);
plot(real(Q),imag(Q));


clear nG
nG(:,:,1) = real(Q);
nG(:,:,2) = imag(Q);

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




%%%%%%%%%%%%%%%%%%%%%%%%
%}