function [curve] = collectContourInformation(fileName,OP)
    % generate nhood
    [x1 x2] = ndgrid(linspace(0,OP.RAD,OP.nRAD),linspace(OP.THETA(1),OP.THETA(2),OP.nTHETA));
    nhoodSZ = size(x1);
    H = [x1(:).*cos(x2(:)) x1(:).*sin(x2(:))]';

    % read and filter the image
    if ischar(fileName)
        I = double(imread(fileName))/255;
    else
        I = fileName;
    end
    fI = imfilter(I,fspecial('gaussian',[11 11],3));
    
    % get level sets
    [curve] = getLevelContours(fI,OP.levels);
    
    % remove non closed curves
    rmidx = [];
    for e = 1:numel(curve)
        if ~all(curve(e).data(:,1) == curve(e).data(:,end))
            rmidx(e) = 1;
        end
    end
    curve(find(rmidx)) = [];
    
    % remove based on length
    rmidx = [];
    for e = 1:numel(curve)
        if curve(e).length > OP.filter_area(2) | curve(e).length < OP.filter_area(1)
            rmidx(e) = 1;
        end
    end
    curve(find(rmidx)) = [];
    
    % display curves
    if OP.displayCurve
        imshow(I,[])
        hold on
        for e = 1:numel(curve)
            plot(curve(e).data(1,:),curve(e).data(2,:),'b');
            drawnow
        end
    end
    
    % generate the normal and tangent space for the curve
    for e = 1:numel(curve)
        curve(e).E = getNormalsAndTangent(curve(e).data,OP.waveletPara);
    end
    
    % sample the curve(s) along the T-N basis space
    for e = 1:numel(curve)
        curve(e).S = sampleImage(curve(e),H,I,nhoodSZ);
    end
    
    % generate curve segments
    for e = 1:numel(curve)
        curve(e).segs = generateCurveSegments(curve(e).data,OP.segmentLength);
    end
    
    %{
    imshow(I,[]);
    hold on;
    for e = 1:numel(curve)     
        hold on
        plot(curve(e).data(1,:),curve(e).data(2,:))
        SAM =10;
        for e1 = 1:SAM:size(curve(e).data,2)
            quiver(curve(e).data(1,e1),curve(e).data(2,e1),curve(e).E(e1,1,1),curve(e).E(e1,1,2),'r');
            quiver(curve(e).data(1,e1),curve(e).data(2,e1),curve(e).E(e1,2,1),curve(e).E(e1,2,2),'g');
            drawnow
        end
    end
    %}
    
    
    %{
    for e = 1:numel(curve)
        for e1 = 1:size(curve(e).segs,1)
            imshow(I,[]);
            hold on
            plot(curve(e).data(1,:),curve(e).data(2,:),'b');
            plot(curve(e).segs(e1,:,1),curve(e).segs(e1,:,2),'r');
            drawnow
            hold off
        end
    end
    %}
end



% generate the normal and tangent space
function [E] = getNormalsAndTangent(segment,S)
    sz = size(segment);
    J = [segment';segment';segment'];
    % calculate curvature
    d1X1 = cwt(J(:,1),S,'gaus1');
    d1X2 = cwt(J(:,2),S,'gaus1');            
    T = cat(3,d1X1,d1X2);
    L = sum(T.*T,3).^.5;
    T = bsxfun(@times,T,L.^-1);
    N = cat(3,T(:,:,2),-T(:,:,1));
    N = squeeze(N)';
    N = N(:,sz(2)+1:sz(2)+sz(2));
    T = squeeze(T)';
    T = T(:,sz(2)+1:sz(2)+sz(2));
    E = cat(3,permute(T,[2 1]),permute(N,[2 1]));
    E = permute(E,[2 3 1]);
end

% sample the image
function [S] = sampleImage(curve,nhood,image,sz)
    % for each point
    for e = 1:size(curve.data,2)
        tmpH = bsxfun(@plus,squeeze(curve.E(:,:,e))*nhood,curve.data(:,e));
        tmpS = ba_interp2(image,tmpH(1,:),tmpH(2,:));
        S(:,:,e) = reshape(tmpS,sz);
        %{
        imshow(image,[]);
        hold on
        plot(tmpH(1,:),tmpH(2,:),'.');
        imshow(tmpS,[]);
        drawnow
        hold off
        %}
    end    
end

% generate the normal and tangent space
function [segs] = generateCurveSegments(segment,S)
    sz = size(segment);
    J = [segment';segment';segment'];
    tmp1 = im2col(J(:,1),[S 1]);
    tmp2 = im2col(J(:,2),[S 1]);
    sz = size(segment);
    tmp1 = tmp1(:,sz(2)+1:sz(2)+sz(2));
    tmp2 = tmp2(:,sz(2)+1:sz(2)+sz(2));
    segs = cat(3,tmp1,tmp2);
    segs = permute(segs,[2 1 3]);
    segs = permute(segs,[2 3 1]);
end


%{
    fileName = SET{1}{1};
    OP.levels = linspace(0,1,10);
    OP.filter_area = [500 3000];
    OP.displayCurve = 0;
    OP.waveletPara = 5;
    OP.RAD = 30;
    OP.nRAD = 30;
    OP.THETA = [-pi pi];
    OP.nTHETA = 200;
    OP.segmentLength = 31;
    curve = collectContourInformation(fileName,OP);
%}
