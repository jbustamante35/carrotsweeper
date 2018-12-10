function [nI] = newRecifiedImage(I,tform,resEstimate)

    if ischar(I)
        I = double(imread(I))/255;
    end

    
    
     for r = 1:4
        I(1,:,1) = 1;
        I(1,:,2) = 0;
        I(1,:,3) = 0;
        I = imrotate(I,90);
    end
    
    
    if ~isempty(tform)

        % approx number of cm for the image
        W = round(.5*(size(I,1).*resEstimate^-1));
        % get number of points 
        NIP = round(2*W*resEstimate);
        % get space in real world co-dinates
        [w1,w2] = ndgrid(linspace(-W,W,NIP),linspace(-W,W,NIP));
        % transform to real world
        XI = tform.transformPointsForward([w2(:),w1(:)]);
        % interpolate real world
        nI = [];
        for k = 1:size(I,3)
            nI(:,k) = ba_interp2(double(I(:,:,k)),XI(:,1),XI(:,2));
        end
        % reshape to color
        nI = reshape(nI,[NIP NIP 3]);
        
        % normalize to [0 1];
        %nI = nI / 255;

        %{
        for k = 1:3
            tmp = nI(:,:,k);
            m = min(tmp(:));
            tmp = tmp - m;
            M = max(tmp(:));
            tmp = tmp / M;
            nI(:,:,k) = tmp;
        end
        %}
    else
        nI = I;
    end
    
end