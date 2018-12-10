function [M] = contourRefine(I,areaTHRESH)
    I = imresize(I,.5);
    %%%%%%%%%%    
    %%% normalize the image
    I = normalize(I);
    
    %%%%%%%%%%
    % obtain the threshold value        
    thresh = graythresh(I);
    %%%%%%%%%%
    % obtain the gradient
    [d1 d2] = gradient(I);
    d = d1.^2 + d2.^2;
    % obtain the entropy
    en = entropyfilt(I,ones(7));
    en = en / max(en(:));
    % obtaint the gradient of the entropy
    [e1 e2] = gradient(en);
    ee = e1.^2 + e2.^2;
    % obtain the edge of image
    d = edge(imfilter(I,fspecial('gaussian',[5 5],7)));
    % obtain the edge of entropy
    en = edge(en);
    % obtain the OR of image edge and entropy edge
    EG = en + d;
    %%%%%%%%%%%%%%%%%%%%
    % loop over parameters
    N = 10;
    NN = 5;
    eN = 4;
    percent = N/255;
    threshVALUES = linspace(thresh-percent,thresh+percent,NN);
    erodeVALUES = -eN:eN;
    for e = 1:numel(threshVALUES)
        %%%%%%%%%%%%%%%%%%%%
        % threshol
        M = I < threshVALUES(e);
        %%%%%%%%%%%%%%%%%%%%
        % loop over erode values
        for i = 1:numel(erodeVALUES)
            %%%%%%%%%%%%%%%%%%%%
            % create structure element
            ST = strel('disk',abs(erodeVALUES(i)));
            if erodeVALUES(i) < 0
                Mn = imerode(M,ST);
            else
                Mn = imdilate(M,ST);
            end
            %%%%%%%%%%%%%%%%%%%%
            % obtain the gradient along the edge
            E = edge(Mn);
            idx = find(E);
            % integrate edge along contour(s)
            eg = sum(d(idx));
            % integrate gradient of entropy along contour(s)
            H = sum(ee(idx));
            % integrate the image along the contour
            Ii = sum(I(idx));
            % integrate along the OR edge
            Gi = sum(EG(idx));
            V(e,i) = Gi*Ii^-1;
        end
    end
    
    [Vmax midx] = max(V(:));
    [ridx cidx] = ind2sub(size(V),midx);
    
    %%%%%%%%%%%%%%%%%%%%
    M = I < threshVALUES(ridx);
    ST = strel('disk',abs(erodeVALUES(cidx)));
    if erodeVALUES(i) < 0
        M = imerode(M,ST);
    else
        M = imdilate(M,ST);
    end
    
    %%%%%%%%%%
    % auto fill holes    
    M = imfill(M,'holes');
    %%%%%%%%%%
    %%% remove small objects
    M = bwareaopen(M,areaTHRESH);
    M = imresize(M,2);
end