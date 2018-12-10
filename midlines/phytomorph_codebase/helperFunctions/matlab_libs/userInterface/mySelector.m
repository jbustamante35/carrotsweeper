function [pS selV] = mySelector(I,pS,gP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % will call point selector with image and snap to the nearest point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I           := image
    %           ps          := point set
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           ps          := selected point set
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    warning off;
    nf = figure;
    imshow(I,[]);
    
    hold on;
    plot(pS(:,2),pS(:,1),'r.');
    if ~isempty(gP); plot(gP(:,2),gP(:,1),'g*');end
    [x1 x2 V] = impixel();
    idx = [];
    for p = 1:numel(x1)
        delta = pS - repmat([x2(p) x1(p)],[size(pS,1) 1]);
        delta = sum(delta.*delta,2);
        [u,idx(p)] = min(delta);
    end
    selV = zeros(size(pS,1),1);
    selV(idx) = 1;
    pS = pS(idx,:);
    
    plot(pS(:,2),pS(:,1),'mo');
    pause(1);
    close(nf);
end