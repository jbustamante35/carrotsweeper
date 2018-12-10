function [iten2] = boundaryManifold(ten1,ten2,GENS,POP,disp)
        %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % match ten1 to ten2
        % assume inf points on ten2 manifold
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%
        polyNum = 10;
        polyDeg = 3;        
        %%%
        % obtain sample spacing of 1
        X = 0:size(ten1,1)-1;
        Y = ones(size(ten1,1),1);
        %%%
        % fit the init spline to the spacing 
        splineT = spap2(polyNum,polyDeg,X,Y);        
        %%%
        % remove the coeffs 
        xo = [splineT.coefs];
        %%%
        % create spacing domain for ten2
        dv = diff(ten2,1,1);
        para = [0;cumsum(sum(dv.*dv,2).^.5)];
        %%%
        % set the optimization parameters
        options = psooptimset('PopInitRange',[zeros(1,size(xo,1));ones(1,size(xo,1))],'PopulationSize',POP,'Display','on','Generations',GENS,'Vectorized','off','CognitiveAttraction',1.5,'SocialAttraction',.5);    
        % run the optimization
        xo = pso(@(xo)corrObj(ten1,ten2,para,splineT,xo,disp),size(xo,2),[],[],[],[],[],[],[],options);
        %%%
        % eval the obtained coeffs
        splineT.coefs = xo;
        xs = fnval(splineT,Y);    
        xs = cumsum(xs');
        xs = xs/xs(end);
        xs = xs*para(end);
        iten2 = interp1(para,ten2,xs);
end