function [K] = geoKur(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute geodesic curvature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script         
    %                   := para.sig         -> sigma for gaussian filter
    %                   := para.gradPara    -> gradient configure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           cM      := corner map       -> corner strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % init OUT vars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        K = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S = para.scales.value;
        imS = para.resize.value;
        sz = size(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resize image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (imS~=1)
        % size image
        Iorg = I;
        I = imresize(I,imS);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first order terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [dFd1 dFd2] = gradient(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % second order terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [d1d1 d1d2] = gradient(dFd1);
        [d2d1 d2d2] = gradient(dFd2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each scale
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:numel(S,2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % construct filter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % disk size
            PDSZ = 2*S(k) + S(k);
            % construct fiter            
            h1 = fspecial('disk', PDSZ);
            h2 = fspecial('gaussian',size(h1),S(k));
            h = h1.*h2;
            % filter first derivative
            tdFd1 = imfilter(dFd1,h,'symmetric');
            tdFd2 = imfilter(dFd2,h,'symmetric');
            % filter second derivative
            td1d1 = imfilter(d1d1,h,'symmetric');
            td1d2 = imfilter(d1d2,h,'symmetric');
            td2d1 = imfilter(d2d1,h,'symmetric');
            td2d2 = imfilter(d2d2,h,'symmetric');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % curve co-ordinate system
            curNOR = cat(3,tdFd1,tdFd2);
            curTAN = cat(3,-tdFd2,tdFd1);
            lam = coVar(curTAN,curTAN,curNOR);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % resize rsults
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (imS~=1)
                lam = imresize(lam,sz);
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store curvature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            K = cat(3,K,lam);
    end


end


function [v] = coVar(X,Y,Z)
    %dX/dY along Z
    cvd = [];
    for i = 1:size(X,3)
        [d1 d2] = gradient(X(:,:,i));
        cvd = cat(3,cvd,Y(:,:,1).*d1 + Y(:,:,2).*d2);
    end
    v = dot(cvd,Z,3);
end





