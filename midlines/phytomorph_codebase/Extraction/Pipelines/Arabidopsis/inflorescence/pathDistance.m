function [d] = pathDistance(p1,p2)


    

    NP = 300;
    d1 = cumsum([0;sum(diff(p1,1,1).^2,2).^.5]);
    s1 = interp1(d1,p1,linspace(0,d1(end),300));
    
    NP = 300;
    d2 = cumsum([0;sum(diff(p2,1,1).^2,2).^.5]);
    s2 = interp1(d2,p2,linspace(0,d2(end),300));
    
    
 
    delta = s2 - s1;
    d(1) = mean(sum(delta.*delta,2).^.5,1);
    d(3) = std(sum(delta.*delta,2).^.5,1,1);
    
    d(2) = abs(d1(end) - d2(end));
    
    
    d = prod(d);
    
    
    
       
    
    
    plot(s1(:,2),s1(:,1),'k')
    hold on
    plot(s2(:,2),s2(:,1),'m')
    title(num2str(d));
    
    waitforbuttonpress
    
end