function [v] = flf(x,p)
    if size(p,1) > numel(x)
        x = x*ones(size(p,1),1);
    end
    for e = 1:size(p,1)
         v(e,:) = p(e,1).*(1 + exp(-p(e,2)*(x(e,:) - p(e,3)))).^-p(e,4);
    end
end


%{
    p = [2 .02 400 1.01];
    x = linspace(0,1000,1000);
    v = flf(x,p);
    plot(x,v);



    [X,T] = position_flf(.1,p,1,.9);
    plot(T,X);

    p = [[2 .02 400 1.01];[2.5 .02 400 1.01]];
    x = linspace(0,1000,1000);
    v = flf(x,p);
    plot(x,v);

    [X,T] = position_flf([.1;.1],p,1,.9);
    plot(T,X);



  %}