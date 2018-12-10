function [curve] = trackFromPointAtGradient_carrot(I,P,initD,maxStep,RHO,RAD,pointDensity,wsigma)
    x = P(1);
    y = P(2);
    %pointDensity = [30 300];
    %RHO = 30;
    %RAD = pi;
    
    %(B,curve(e).data(:,tipIDX{e}),t1,t2,5000,20,pi/2,[20 200],.7);
    
    %%%%%%%%%%%%%%%%%%
    % the width of the belief or the momentum is a function of percent 
    % of cutoff
    %wsigma = .3;
    maxWidth = 60*pi/180;
    width = @(delta)maxWidth*normpdf(delta,0,wsigma)/normpdf(0,0,wsigma);
    %%%%%%%%%%%%%%%%%%
    % function for radial belief
    k = 1.5;
    alpha = 1;
    scale = 10;
    radial = @(x)wblpdf(x*scale^-1,alpha,k);
    %%%%%%%%%%%%%%%%%%
    % function for angular belief - 
    angle = @(x,delta)normpdf(x,0,width(delta));
    %%%%%%%%%%%%%%%%%%
    % belief function
    func = @(rad,rho,delta)angle(rad,delta).*radial(rho);

    
    %{
    t1 = ba_interp2(g1,x,y);
    t2 = ba_interp2(g2,x,y);
    N = [t1 t2];
    N = N / norm(N);
    T = [-N(2) N(1)];
    initD = [T;N];
    %}
    
    
    T = goT();
    T.setWfunction(func);
    T.setNhoodRho(RHO);
    T.setNhoodRad(RAD);
    T.setNhoodDensity(pointDensity);
    T.generateH();
    T.setPosition([x;y]);
    T.setImage(I);
    T.setDirection(initD);
    stopFunction = @(path) size(path,2) < maxStep & path(1,end) > 300;
    %T.walk(maxStep);
    T.walkUntil(stopFunction);
    T.reparameterizeCurve();
    curve = T.position;

end