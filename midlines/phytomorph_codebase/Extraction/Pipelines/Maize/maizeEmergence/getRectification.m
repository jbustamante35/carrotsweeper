function [rec] = getRectification(imageFile,autoFLAG,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rectify the image based on the checker board
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting image rectification \n']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image and make gray scale
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image file
    I = imread(imageFile);
    % convert RGB to Gray
    G = rgb2gray(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the corners
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CM = cornermetric(G,'Harris','SensitivityFactor',0.01);
    % if not autoFLAG then crop by hand
    if ~autoFLAG
        % make boundingbox via hand
        figure;
        imshow(G);
        h = impoly;
        position = wait(h);
        CBM = poly2mask(position(:,1),position(:,2),size(G,1),size(G,2));
        %CBM = imerode(CBM,strel('square',50));
    else
        paperMASK = getCheckerBoardMask(I);
        CBM = imfill(paperMASK,'holes') - paperMASK;
        CBM = imerode(CBM,strel('disk',50,0));
    end
    
    
    % make the corner mask
    %CM = bindVec(CM(:));
    CMM = CM >0.000005;
    %CMM = CM > .95;
    CMM = (CM == imdilate(CM,strel('disk',100))) & (CM > 0.000001);
    CMM = logical(CMM.*CBM);
    % find the centroids of the corners
    R = regionprops(CMM,'Centroid');
    sP = [];
    for e = 1:numel(R)
        sP = [sP;R(e).Centroid];
    end
    if disp
        close all;
        imshow(I,[]);hold on;
        for e = 1:numel(R)
            plot(R(e).Centroid(1),R(e).Centroid(2),'r.');
        end
    end
    usP = mean(sP,1);

    
    
    rsP = [sP(:,2) -sP(:,1)];
    [J pidx(1)] = max(sum(sP.*sP,2));
    [J pidx(2)] = min(sum(sP.*sP,2));


    tanV = sP(pidx(1),:) - sP(pidx(2),:);
    tanV = tanV /norm(tanV);
    norV = [-tanV(2) tanV(1)];

    osP = bsxfun(@minus,sP,usP);
    alongVec = osP*norV';
    [alongVec sidx] = sort(alongVec);
    pidx(3) = sidx(1);
    pidx(4) = sidx(end);

    if disp
        quiver(usP(1),usP(2),tanV(1),tanV(2),1000);
        quiver(usP(1),usP(2),-tanV(1),-tanV(2),1000);
        quiver(usP(1),usP(2),norV(1),norV(2),1000);
        quiver(usP(1),usP(2),-norV(1),-norV(2),1000);
    end

    schIdx = [2 1 3 4];
    pidx = pidx(schIdx);

    if disp
        plot(usP(1),usP(2),'g.');
        CL = {'ro' 'bo' 'go' 'ko'};
        for e = 1:numel(pidx)
            plot(sP(pidx(e),1),sP(pidx(e),2),CL{e});
        end
    end

    sourcePoints = sP(pidx,:);
    targetPoints = sP(pidx,:);
    targetPoints(1,1) = sourcePoints(4,1);
    LEN = targetPoints(4,2) - targetPoints(1,2);
    targetPoints(3,:) = [targetPoints(4,1) + LEN targetPoints(4,2)];
    targetPoints(2,:) = [targetPoints(3,1) targetPoints(1,2)];

    tmp = sourcePoints(2,:);
    sourcePoints(2,:) = sourcePoints(3,:);
    sourcePoints(3,:) = tmp;
    if disp
        for e = 1:size(sourcePoints,1)
            plot(sourcePoints(e,1),sourcePoints(e,2),'r*')
            text(sourcePoints(e,1),sourcePoints(e,2),['s' num2str(e)]);
            plot(targetPoints(e,1),targetPoints(e,2),'k*')
            text(targetPoints(e,1),targetPoints(e,2),['t' num2str(e)]);
        end
    end

    rec = cp2tform((sourcePoints),(targetPoints), 'projective');
    if disp
        G = imtransform(G,rec);
        figure
        imshow(G,[]);
    end
    
    fprintf(['ending image rectification \n']);
end