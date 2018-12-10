function [curve] = alignCurves(curve)
    curve = interpolateCurves(curve,1000);
    for i = 1:numel(curve)
        
    end
end

function [d] = curveDistance(X,Y)
    
    [out] = cwtK_closed(J,para)
end

function [curve] = interpolateCurves(curve,N)
    for e = 1:numel(curve)
        d = diff(curve(e).data,1,2);
        l = sum(d.*d,1).^.5;
        l = cumsum([0 l]);
        curve(e).data = interp1(l,curve(e).data',linspace(0,l(end),N));
    end
end


