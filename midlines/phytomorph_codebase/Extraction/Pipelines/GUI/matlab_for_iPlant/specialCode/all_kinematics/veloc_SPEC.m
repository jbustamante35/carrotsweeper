function vel = veloc_SPEC(p,x,fixed)
    try
        if nargin == 2
            fixed.k = ones(size(p));
            fixed.v = zeros(size(p));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vo,vf,n,k,xo
        % p(1) = value of velocity @ xo=p(5)
        % p(2) = growth rate of the root
        % p(3) = shape parameter
        % p(4) = shape parameter
        % p(5) = position for p(1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        p = p.*fixed.k + fixed.v;

        vel = p(1)*p(2).*(p(1)^p(3) + (p(2)^p(3) - p(1)^p(3))*exp(-p(4)*(x - p(5)))).^-(1/p(3));    

        if any(isnan(vel(:)))
            vel = 100000*ones(size(vel));
        end
    catch ME
        stop= 0;
    end
end