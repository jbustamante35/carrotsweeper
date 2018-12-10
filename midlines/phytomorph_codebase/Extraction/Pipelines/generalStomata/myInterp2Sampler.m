function [sI] = myInterp2Sampler(I,indexPosition,Domain,rsz)
    dispDebug = false;
    if numel(indexPosition) == 1
        % generate n-sub index
        [xPos(2,:),xPos(1,:)] = ind2sub(size(I),indexPosition);
    else
        xPos = indexPosition;
    end
    % displace the domain
    displaceDomain = bsxfun(@plus,Domain,xPos');
    % create the storage vector
    sI = zeros(size(displaceDomain,1),size(I,3));
    % loop over the color space
    for e = 1:size(I,3)
        sI(:,e) = ba_interp2(I(:,:,e),displaceDomain(:,1),displaceDomain(:,2),'cubic');
        %sI(:,e) = zscore(sI(:,e),1,1);
        %sI(:,e) = bsxfun(@minus,sI(:,e),mean(sI(:,e),1));
    end
    % resize if needed
    if nargin == 4
        sI = reshape(sI,[rsz size(sI,2)]);
    end
    if dispDebug
        imshow(I,[]);
        hold on
        plot(displaceDomain(:,1),displaceDomain(:,2),'g.');
        imshow(sI,[])
        
        hold off
        drawnowxPos
        
        
        
        
    end
end