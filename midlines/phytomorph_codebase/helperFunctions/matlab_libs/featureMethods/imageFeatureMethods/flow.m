function [P] = flow(V,P,n,s,idisp)
    dispV = 1;
    for i = 1:n
        
        vec1 = ba_interp2(V(:,:,1),P(end,1),P(end,2));
        vec2 = ba_interp2(V(:,:,2),P(end,1),P(end,2));
        
        vec = [vec1 vec2];
        
        vec = vec * norm(vec)^-1;
        vec = vec * s;
        P = [P;P(end,:) + vec];
        
        if dispV;imshow(idisp,[]);hold on;plot(P(:,1),P(:,2),'r');hold off;drawnow;end
    end
end