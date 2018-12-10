function [D] = midObjFunc(set1,set2,para,disp,spline,xo)
    % match set1 to set2
    % spline the ds along the set2
    % interpolate the ds along set2 at the number of points from set1
    % 
    spline.coefs = xo';
    xo = fnval(spline,0:size(set1,1)-1);
    
    O = xo;
    xo = cumsum(xo);
    xo = xo/xo(end);
    xo = xo*para(end);
    set2i = interp1(para,set2,xo);
    midP = (set2i - set1)/2 + set1;
    
    connVec = (set2i - set1);
    
    T = diff(midP,1,1);
    
    connVec(end,:) = [];
    
    for i = 1:size(T,1)
        T(i,:) = T(i,:)/norm(T(i,:));
        connVec(i,:) = connVec(i,:)/norm(connVec(i,:));    
    end
    
    D = sum(abs(sum(T.*connVec,2)));
    
    
    if disp
        SKIP = 10;
        plot(set2(:,2),-set2(:,1),'g')
        hold on
        plot(set1(:,2),-set1(:,1),'r')
        for i = 1:SKIP:size(set1,1)
            plot([set2i(i,2) set1(i,2)],-[set2i(i,1) set1(i,1)])    
        end
        quiver(midP(1:SKIP:end-1,2),-midP(1:SKIP:end-1,1),connVec(1:SKIP:end,2),-connVec(1:SKIP:end,1),'r')
        quiver(midP(1:SKIP:end-1,2),-midP(1:SKIP:end-1,1),T(1:SKIP:end,2),-T(1:SKIP:end,1),'c')
        plot(midP(1:SKIP:end,2),-midP(1:SKIP:end,1),'g')
        drawnow
        hold off
    end
end

