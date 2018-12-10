function [] = genDomain(L,nL,H,nH,S,nS,m,nE,disp)
    M = (L - 2*S)/2;
    TH = linspace(0,pi,nE);
    delta = L/2;
    P1 = [linspace(0,L,nL);linspace(H,H,nL)];
    P2 = [linspace(L,L,nH);linspace(H,0,nH)];
    P3 = [linspace(L,L-S,nS);linspace(0,0,nS)];
    P4 = [M*cos(TH)+delta;m*sin(TH)];
    P5 = [linspace(S,0,nS);linspace(0,0,nS)];
    P6 = [linspace(0,0,nH);linspace(0,H,nH)];
    dB = [P1 P2 P3 P4 P5 P6];
    if disp
        plot(dB(1,:),dB(2,:),'r');
    end
end

%{  
    genDomain(100,100,10,10,5,20,3,100,1);
%}