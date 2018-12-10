function [stop] = myPF1(data,x,state)
    nf = numel(x.bestx)/2;
    [yv,xv] = hist([data-1;data;data+1],500);
    plot(xv,yv/sum(yv))
    XV = linspace(-1,2,500);
    hold on
    for k = 1:nf
        Y =  normpdf(XV,x.bestx(k),x.bestx(k+nf));
        Y = Y / sum(Y);
        plot(XV,Y)
        hold on
    end
    
    hold off
    drawnow
    stop = false;
end