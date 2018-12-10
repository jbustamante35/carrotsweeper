function [curve] = collectContourInformation_ver2(fileName,OP)
    % generate save information
    [pth nm ext] = fileparts(fileName);
    
    % generate nhood
    [x1 x2] = ndgrid(linspace(0,OP.RAD,OP.nRAD),linspace(OP.THETA(1),OP.THETA(2),OP.nTHETA));
    nhoodSZ = size(x1);
    H = [x1(:).*cos(x2(:)) x1(:).*sin(x2(:))]';

    % read and filter the image
    I = double(imread(fileName))/255;
    fI = imfilter(I,fspecial('gaussian',[11 11],2),'replicate');
    fI = bindVec(fI);
    
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
    
    
    CM = generateContainmentMap(curve);
    rmidx = ~all(CM==1,2);
    curve(rmidx) = [];
    
    curve = arcLengthNormalizeCurves(curve);
    
    for e = 1:numel(curve)
        % set the tmp curve
        tmpCurve = curve(e);
        % generate the normal and tangent space for the curve
        tmpCurve.E = getNormalsAndTangent(tmpCurve.data,OP.waveletPara);
        % sample the curve(s) along the T-N basis space
        tmpCurve.S = sampleImage(tmpCurve,H,I,nhoodSZ);
        % generate curve segments
        tmpCurve.segs = generateCurveSegments(tmpCurve,OP.segmentLength);
        
        %tmpCurve = decomposeData(tmpCurve);
        
        % label from which image the curve comes from
        tmpCurve.imageName = fileName;
        % save the data
        oFileName =[OP.savePath strrep(pth,'/','_') nm '_curve_' num2str(e) '.mat'];
        save(oFileName,'tmpCurve');
        
        
        
        % display curves
        if OP.displayCurve
            close all
            imshow(fI,[])
            hold on
            for e = 1:numel(curve)
                plot(tmpCurve.data(1,:),tmpCurve.data(2,:),'b');
                quiver(tmpCurve.data(1,:)',tmpCurve.data(2,:)',squeeze(tmpCurve(e).E(1,2,:)),squeeze(tmpCurve(e).E(2,2,:)));
                drawnow
            end
        end  

        
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
function [segs] = generateCurveSegments(curve,segmentSize)
    halfBlock = (segmentSize-1)/2;
    segment = curve.data;
    sz = size(segment);
    J = [segment';segment';segment'];
    tmp1 = im2col(J(:,1),[segmentSize 1]);
    tmp2 = im2col(J(:,2),[segmentSize 1]);    
    tmp1 = tmp1(:,sz(2)+1-halfBlock:sz(2)+sz(2)-halfBlock);
    tmp2 = tmp2(:,sz(2)+1-halfBlock:sz(2)+sz(2)-halfBlock);
    segs = cat(3,tmp1,tmp2);
    segs = permute(segs,[2 1 3]);
    segs = permute(segs,[2 3 1]);
    for e = 1:size(segs,3)
        segs(:,:,e) = bsxfun(@minus,segs(:,:,e),mean(segs(:,:,e),1));
        segs(:,:,e) = (curve.E(:,:,e)'*segs(:,:,e)')';
    end
    %{
    for e = 1:size(segs,3)
        plot(segs(:,1,e),segs(:,2,e))
        hold on
    end
    %}
end

% PCA on curve data
function [curve] = decomposeData(curve)
    % gather stats on the curves points
    [curveSim curveC curveU curveE curveL ERR curveLAM] = PCA_FIT_FULL(curve.data',2);
    % gather stats on the curves patches
    tmpData = reshape(curve.S,[size(curve.S,1)*size(curve.S,2) size(curve.S,3)]);
    [patchSim patchC patchU patchE patchL ERR patchLAM] = PCA_FIT_FULL(tmpData',5);
    % gather stats on the curves segs
    tmpData = reshape(curve.segs,[size(curve.segs,1)*size(curve.segs,2) size(curve.segs,3)]);
    [segsSim segsC segsU segsE segsL ERR segsLAM] = PCA_FIT_FULL(tmpData',5);
    
end



function [curves] = arcLengthNormalizeCurves(curves)
    for e = 1:numel(curves)
        dL = diff(curves(e).data,1,2);
        dL = sum(dL.*dL,1).^.5;
        L = cumsum([0 dL]);
        curves(e).data = interp1(L',curves(e).data',linspace(0,L(end),round(L(end)))')';
    end
end


function [CM] = generateContainmentMap(tmpC)
    CM = [];
    for i = 1:numel(tmpC)        
        PTL = [];        
        for j = 1:numel(tmpC)        
            PTL = [PTL;tmpC(j).data(:,1)'];         
        end
        CM(i,:) = inpoly(PTL,tmpC(i).data(:,1:10:end)');        
    end
end
%{
    fileName = SET{1}{1};
    fileName = outFile;
    OP.levels = linspace(0,1,100);
    OP.filter_area = [500 3000];
    OP.displayCurve = 0;
    OP.waveletPara = 5;
    OP.RAD = 30;
    OP.nRAD = 30;
    OP.THETA = [-pi pi];
    OP.nTHETA = 200;
    OP.segmentLength = 31;
    OP.savePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/development/maizeContour/';
    parfor s = 1:numel(SET)
        fileName = SET{s}{1};        
        curve = collectContourInformation_ver2(fileName,OP);
    end

%}
