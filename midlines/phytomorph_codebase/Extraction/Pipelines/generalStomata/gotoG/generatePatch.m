function [] = generatePatch()


%%
PZ = [100 100];
[n1,n2] = ndgrid(-PZ(1):PZ(1),-PZ(2):PZ(2));
rho = (n1.^2 + n2.^2).^.5;
theta = atan2(n2,n1);

k1 = 5;
phi1 = pi/5;
v1 = 1;
sig1 = 0;


omega = @(k,phi,v,sig,theta)k*(phi + theta + v*sig);

F = exp(i*omega(k1,phi1,v1,sig1,theta));
close all
imshow(real(F),[]);

%% sweep phi
close all
phiS = linspace(-pi,pi,100);
for e = 1:numel(phiS)
    F = exp(i*omega(k1,phiS(e),v1,sig1,theta));
    imshow(F,[]);
    drawnow
    pause(.2)
end
%% sweep sig
close all
sigS = linspace(-pi,pi,100);
for e = 1:numel(phiS)
    F = exp(i*omega(k1,phi1,v1,sigS(e),theta));
    imshow(F,[]);
    drawnow
    pause(.2)
end
%%

k_vec = [2 5 4];
phi_vec = [0 pi/4 -3*pi/4];
v_vec = [1 1 1];
sig_vec = [0 0 0];


[wave] = buildWave(k_vec,v_vec,sig_vec,phi_vec,theta);
close all
imshow(real(wave),[])
%%
rawI = imread(sorFileList{4});
toOp = rawI;

R = [0 40];
N = [41 250];
NF = 10;
P = 0;
mag = 1;
func = [];
disp = 0;
[retVec1,retVec2,rawI] = gogo_bugEye(rawI,toOp,R,N,NF,mag,P,disp,func);
%%
PZ = [R(2) R(2)];
[n1,n2] = ndgrid(-PZ(1):PZ(1),-PZ(2):PZ(2));
rho = (n1.^2 + n2.^2).^.5;
theta = atan2(n1,n2);
wave_tensor = squeeze(retVec1);
v_vec = ones(size(wave_tensor));
sig_vec = pi/2*ones(size(wave_tensor));
[wave] = buildWave(wave_tensor,v_vec,sig_vec,theta,rho,R,NF,N);
end