function [DIS] = minPhaseDifference(cPhi,k1,phi1,k2,phi2,disp)
    %{
    phi2 = phi2 * k2.^-1;
    phi1 = phi1 * k1.^-1;

    dPhi = (phi2 - phi1);
    expectedK = .5*(k1^-1 - k2^-1);
    DIS = (2*pi).^2 * expectedK + dPhi;
    %}
    
    
    % OLD FOR ONE LOCATION
    
    phi2 = phi2 * k2 ^-1;
    phi1 = phi1 * k1 ^-1;
    
    
    
    %{
    deltaTheta = (phi2 - phi1);
    choiceTheta = -sign(deltaTheta).*2*pi + sign(deltaTheta).*deltaTheta;
    choiceVector = [deltaTheta;choiceTheta];
    [~,idx] = min(abs(choiceVector),[],1);
    DIS = choiceVector(idx);
    %}
    

    if nargin == 5
        disp = false;
    end
    mag = 1;
    n1 = (-mag*k1:mag*k1) - 1;
    n2 = (-mag*k2:mag*k2) - 1;
    
    
    %{
    xv = @(n,phiC,k,phi)(n*2*pi + (phiC - phi))*k^-1;
    x1 = xv(n1,cPhi,k1,phi1);
    x2 = xv(n2,cPhi,k2,phi2);
    %}
    
    x1 = (n1*2*pi + cPhi - phi1*k1)*k1^-1;
    x2 = (n2*2*pi + cPhi - phi2*k2)*k2^-1;
    

    nx1 = (k1^-1)*2*pi*n1 - phi1;


    %{
    x1 = (n1*2*pi + cPhi)*k1^-1 - phi1;
    x2 = (n2*2*pi + cPhi)*k2^-1 - phi2;
    %}

    %{
    x1(x1 < -pi | x1 > pi) = [];
    x2(x2 < -pi | x2 > pi) = [];
    %}
    
    %tD = pdist2(x1',x2');
    
    
    tD = zeros(numel(x1),numel(x2));
    for e1 = 1:numel(x1)
        for e2 = 1:numel(x2)
            tD(e1,e2) = (x1(e1) - x2(e2));
        end
    end
    
    [DIS,dIDX] = min(abs(tD(:)));
    DIS = tD(dIDX);
    
    




    if disp
        %{
        phi2 = phi2 * k2 ^-1;
        phi1 = phi1 * k1 ^-1;
%}
        %{
        x1(x1 < -pi | x1 > pi) = [];
        x2(x2 < -pi | x2 > pi) = [];
        %}
        
        TOTN = 1;
        numP = TOTN*1000;
        
        wave1 = @(x)exp(1i*(k1*(x+phi1)));
        wave2 = @(x)exp(1i*(k2*(x+phi2)));
        
        x = linspace(TOTN*-pi,pi,numP);
        
        
        superWave = wave1(x) + wave2(x);
        plot(x,real(wave1(x)),'r--')
        hold on
        plot(x,real(wave2(x)),'b')
        plot(x,superWave,'g')
        plot(x1,wave1(x1),'r*');
        plot(x2,wave2(x2),'b*');
        
        cV = exp(1i*cPhi);
        plot(x,real(cV)*ones(size(x)),'k')

        
        
        title(num2str(DIS));
        [i1,i2] = ind2sub(size(tD),dIDX);
        mX1 = [x1(i1) wave1(x1(i1))];
        mX2 = [x2(i2) wave2(x2(i2))];
        SEG = [mX1;mX2];

        plot(SEG(:,1),real(SEG(:,2)),'k','LineWidth',2);
        axis([-pi pi -2 2]);
        drawnow
        hold off
    end
    
end

%{
    close all
    k1 = 3;
    k2 = 4;
    phi1 = 4*pi/7;
    phi2 = 5*pi/32;
    dP = linspace(-2*pi,2*pi,300);
    for e = 1:numel(dP)
        [DIS(e)] = minPhaseDifference(0,k1,phi1+dP(e),k2,phi2+dP(e),true);
    end


    TH = linspace(-pi,pi,1000);
    w = @(x,k,p)exp(i*(k*(x+p)));
    ps = linspace(-pi/2,pi/2,50);
    k1 = 3;
    k2 = 4;
    p1 = -pi/8;
    p2 = pi/9;
    vel = 1;
    ps1 = p1 + vel*linspace(-2*pi,2*pi,500);
    ps2 = p2 + vel*linspace(-2*pi,2*pi,500);
delta = [];
    for e = 1:numel(ps1)
        y1 = w(TH,k1,ps1(e));
        y2 = w(TH,k2,ps2(e));

        f1 = fft(y1);
        a1 = angle(f1(k1+1));
       
        f2 = fft(y2);
        a2 = angle(f2(k2+1));
    
 A1(e) = a1/k1;
        A2(e) = a2/k2;

        delta(e) = minPhaseDifference(0,k1,a1,k2,a2,true);
        %{
        
        
        plot(TH,y1,'b');
        hold on
        plot(TH,y2,'r');
        drawnow
        pause(.2)
        hold off
        %}
    end





%}