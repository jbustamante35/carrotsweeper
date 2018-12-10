function [  ] = plotTassel( tBin, tBinProps, splines, endpoints, base, circle, lowBranch, maxSpline, topBranch, spikeTip )
%PLOTTASSEL Creates an image of the binary tassel with some features drawn
%on.

    
    colormap('gray')
    image(tBin, 'CDataMapping', 'scaled'); hold on;
    plot(tBinProps.ConvexHull(:,1), tBinProps.ConvexHull(:,2), 'c-')   
    scatter(endpoints(:,1), endpoints(:,2), 100, 'r', '+');

    for spline=1:size(splines, 2);
        S = splines(spline);
        fnplt(S, 2, 'g');
    end
    
    if(exist('maxSpline', 'var'))
        fnplt(maxSpline, 'b');
        coords = ppval(maxSpline, [1 size(maxSpline.coefs, 1)/2 ]);
        plot(coords(1,:), coords(2,:), 'red', 'LineWidth', .75);
    end
    
    if(exist('spikeTip', 'var'))
        fnplt(spikeTip, 'm');
    end
    
    if(exist('topBranch', 'var'))
        scatter(topBranch(1), topBranch(2), 75, 'blue', 'filled')
    end
    
    if ~isnan(circle)
        for c=1:size(circle, 1);
            plot(circle(c,1), circle(c,2), 'yellow')
        end
        scatter(lowBranch(1), lowBranch(2), 50, 'blue', 'filled')
    end

    scatter(base(1), base(2), 25, 'filled', 'ms');
    
    hold off;
end

