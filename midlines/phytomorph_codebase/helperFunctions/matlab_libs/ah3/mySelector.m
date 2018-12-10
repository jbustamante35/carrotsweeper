function [pS] = mySelector(I,pS)
    nf = figure;
    imshow(I,[]);
    hold on
    plot(pS(:,2),pS(:,1),'r.');
    [x1 x2 V] = impixel();
    for p = 1:numel(x1)
        delta = pS - repmat([x2(p) x1(p)],[size(pS,1) 1]);
        delta = sum(delta.*delta,2);
        [u,idx(p)] = min(delta);
    end
    pS = pS(idx,:);
    plot(pS(:,2),pS(:,1),'go');
    waitforbuttonpress;
    close(nf);
end