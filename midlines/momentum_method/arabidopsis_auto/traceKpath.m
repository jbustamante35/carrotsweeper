function [P dist] = traceKpath(P,W1,W2)
    THRESH = 10;
    MAX = 900;
    e = 1;
    flag = 1;
    while flag~=0 && e <= MAX;
        for p = 1:size(P,1)
            d1 = ba_interp2(W1,P(p,1,e),P(p,2,e));
            d2 = ba_interp2(W2,P(p,1,e),P(p,2,e));
            d = [d1 d2];
            d = d / norm(d);
            d = d * .5;
            
            if e ~= 1
                dV = squeeze(P(p,:,e) - P(p,:,e-1));
                if d*dV'< 0
                    d = -d;                
                end
            end
            
            P(p,1,e+1) = P(p,1,e) + d(1);
            P(p,2,e+1) = P(p,2,e) + d(2);
        end
        [J,dP] = gradient(P(:,:,end));
        dist(e+1) = squeeze(sum(dP(2,:).*dP(2,:),2).^.5);
        if dist(e+1) > THRESH;
            flag=0;
        end
        e = e + 1;
    end
end