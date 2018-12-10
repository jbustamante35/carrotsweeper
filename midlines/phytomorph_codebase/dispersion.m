T1 = 2;
T2 = 3;


k1 = 2*pi/(T1);
k2 = 2*pi/(T2);


phi1 = pi/2;
phi2 = -3*pi/5;

deltaPhi = (phi2 - phi1)/k1;
rho = k2/k1;
N = 2;
M = 1000*N;
TH = linspace(-pi*N,pi*N,M);
Y = exp(-i*(rho*TH + deltaPhi));

wave1 = exp(i*(k1*TH+phi1));
wave2 = exp(i*(k2*TH+phi2));
close all
plot(TH,real(Y));
figure;
plot(TH,real(wave1),'r');
hold all
plot(TH,real(wave2),'b');
plot(TH,real(wave1+wave2),'g')
Tnew = 2*pi/rho
%%
close all
N = 2;
M = 100*N;
TH = linspace(-pi*N,pi*N,M);
kModel = 2;
phiModel = pi/2;
x = TH;
y = exp(i*(kModel*x + phiModel))
plot(x,real(y))
%% test constant
cPhi = pi/4 - 7/8;
testPhi = pi/2;
n = 0:10;
kTest = 3;
xv = n*2*pi/kTest + (cPhi - testPhi)/kTest; 
c = exp(i*(n*2*pi + cPhi));
testWave = @(x)exp(i*(kTest*x + testPhi));
testWave(xv)
%%
cPhi = pi/4 - 7/8;
n = 0:10;
kTest = 3;
xv = n*2*pi/kTest + (cPhi - testPhi)/kTest; 
c = exp(i*(n*2*pi + cPhi));
%%
close all
clear all

k1 = 2;
k2 = 3;
TOTN = 1;
numP = TOTN*1000;
x = linspace(TOTN*-pi,pi,numP);
phi1 = pi/2;
phi2 = -3*pi/5;
n1 = (-k1:k1) - 1;
n2 = (-k2:k2) - 1;



cPhi = pi/4 - 7/8;
cPhi = 0;
cV = exp(i*cPhi);


wave1 = @(x)exp(i*(k1*x+phi1));
wave2 = @(x)exp(i*(k2*x+phi2));
superWave = wave1(x) + wave2(x);
plot(x,real(wave1(x)),'r')
hold on
plot(x,real(wave2(x)),'b')
plot(x,superWave,'g')


plot(x,real(cV)*ones(size(x)),'k')


xv = @(n,phiC,k,phi)(n*2*pi + (phiC - phi))*k^-1;
x1 = xv(n1,cPhi,k1,phi1);
x2 = xv(n2,cPhi,k2,phi2);

x1(x1 < -pi | x1 > pi) = [];
x2(x2 < -pi | x2 > pi) = [];

plot(x1,wave1(x1),'r*');
plot(x2,wave2(x2),'b*');

for e1 = 1:numel(x1)
    for e2 = 1:numel(x2)
        tD(e1,e2) = abs(x1(e1) - x2(e2));
    end
end
[DIS,dIDX] = min(tD(:));
title(num2str(DIS));
[i1,i2] = ind2sub(size(tD),dIDX);
mX1 = [x1(i1) wave1(x1(i1))];
mX2 = [x2(i2) wave2(x2(i2))];
SEG = [mX1;mX2];
plot(SEG(:,1),real(SEG(:,2)),'k','LineWidth',2)








