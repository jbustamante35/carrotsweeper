function [rI] = getRectifiedImage(I,tform,cameraShift)

    fprintf(['Starting image rectification. \n']);
    if nargin <= 2
        cameraShift = [];
    end


    if ischar(I)
        I = imread(I);
        sz = size(I);
        %{
        if ~isempty(cameraShift)
            for k = 1:3
                I(:,:,k) = imwarp(I(:,:,k),cameraShift,'OutputView',imref2d(sz(1:2)));
            end
        end
        %}
    end
    %[R, t] = extrinsics(imagePoints, worldPoints, cameraParameters);
    
    %NP = 5000;
    
    %imageBORDER = [1 1;size(I,1) 1;size(I,1) size(I,2);1 size(I,2);1 1];
    %worldimageBORDER = pointsToWorld(cameraParameters,R,t,[imageBORDER(:,1) imageBORDER(:,2)]);
    
    
    
    %W = 25.4*4;
    W = 400;
    NIP = 5700;
    %W = 200;
    %NIP = 200;
    [w1 w2] = ndgrid(linspace(-W,W,NIP),linspace(-W,W,NIP));
    %UD = [min(worldimageBORDER(2,2)),400];
    %LR = [-min(worldimageBORDER(4,2)),min(worldimageBORDER(3,2))];
    %UD = [337 
    %RATIO = (UD(2)-UD(1))/(LR(2)-LR(1));
    %newR = round(RATIO^-1*NP);
    
    %{
    UD = [-284 256];
    LR = [-128 222];
    %}
    
    
    %COL = linspace(LR(1),LR(2),newR);
    %ROW = linspace(UD(1),UD(2),NP);
    
    
    
    
    %dX1 = linspace(max([min(X(:,1)) 1]),min([max(X(:,1)) size(I,2)]),NP);
  %  dX2 = linspace(max([min(X(:,2)) 1]),min([max(X(:,2)) size(I,1)]),NP);
    
    %[X1 X2] = ndgrid(ROW,COL);
    
    %[X1 X2] = ndgrid(ROW,COL);
%    XI = worldToImage(cameraParameters,R,t,[w2(:) w1(:) zeros(numel(w1),1)],'ApplyDistortion',false);
    XI = tform.transformPointsForward([w2(:),w1(:)]);
    %XI = bsxfun(@plus,XI,imageCenterPoint);
    XI1 = reshape(XI(:,1),[NIP NIP]);
    M1 = XI1 >= 1 & XI1 <= size(I,2);
    XI2 = reshape(XI(:,2),[NIP NIP]);
    M2 = XI2 >= 1 & XI2 <= size(I,1);
    M = M1 .*M2;
    % get bounding box with regionn props
    fprintf('Getting bounding box with regionprops.');
    R = regionprops(logical(M),'BoundingBox');
    %XI = bsxfun(@plus,XI,imageCenterPoint);
               
    
    rI = ba_interp2(double(I),XI(:,1),XI(:,2));
    rI = reshape(squeeze(rI),[NIP NIP 3])/255;
    rI = bsxfun(@times,rI,M);
    rI = imcrop(rI,R(1).BoundingBox);

    if ~isempty(cameraShift)
        sz = size(rI);
        for k = 1:3
            rI(:,:,k) = imwarp(rI(:,:,k),cameraShift,'OutputView',imref2d(sz(1:2)));
        end
    end

    fprintf(['Ending image rectification. \n']);
  
    %{
    imshow(rI,[]);
    drawnow
    %}
    
end