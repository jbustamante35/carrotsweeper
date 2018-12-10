function [curve] = alignCurves(curve)
    curve = interpolateCurves(curve,1000);
    
    viewCurves(curve(10),curve(50))
    for i = 1:numel(curve)
        for j = i+1:numel(curve)
            [D(i,j) idx(i,j)] = alignTo(curve(i),curve(j),100);
        end
    end
end


function [D idx] = alignTo(T,S,para)
    % measure source kurvature
    Sk = cwtK_closed(S.data',{para});
    Sk = Sk.K;
    for i = 1:size(S.data,2)
        
        Tk = cwtK_closed(T.data',{para});
        Tk = Tk.K;
        d = Sk - Tk;
        d = sum(d);
        dist(i) = d;
        
        
        viewCurves(S,T);
        
        S.data = circshift(S.data,[0 -1]);
    end
    [D idx] = min(dist);
end



function [d] = curveDistanceK(Sk,TC,para)    
    Yk = cwtK_closed(Y.data',{para});
    Yk = Yk.K;
    d = Yk - Xk;
    d = sum(d);
end

function [d] = curveDistanceK(X,Y,para)
    Xk = cwtK_closed(X.data',{para});
    Xk = Xk.K;
    Yk = cwtK_closed(Y.data',{para});
    Yk = Yk.K;
    d = Yk - Xk;
    d = sum(d);
end

function [curve] = interpolateCurves(curve,N)
    for e = 1:numel(curve)
        d = diff(curve(e).data,1,2);
        l = sum(d.*d,1).^.5;
        l = cumsum([0 l]);
        curve(e).data = interp1(l,curve(e).data',linspace(0,l(end),N))';
    end
end

function [] = viewCurves(X,Y)
    uX = mean(X.data,2);
    X.data = bsxfun(@minus,X.data,uX);
    uY = mean(Y.data,2);
    Y.data = bsxfun(@minus,Y.data,uY);
    plot(X.data(1,:),X.data(2,:),'r');
    hold on
    plot(Y.data(1,:),Y.data(2,:),'b');
    DF = X.data - Y.data;
    for e = 1:5:size(DF,2)
        plot([Y.data(1,e) Y.data(1,e)+DF(1,e)],[Y.data(2,e) Y.data(2,e)+DF(2,e)],'k')
    end
    drawnow
    hold off
end

