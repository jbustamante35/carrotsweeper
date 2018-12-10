function [] = zoomSparkleLock_QR(e,fileName,boxSequence,zoomSequence,map_Q,beta_Q,WINDOW_Q,SPboxSequence,SPzoomSequence,map_S,beta_S,WINDOW_S,boxSequence_F,zoomSequence_F,map_F,beta_F,WINDOW_F)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = double(imread(fileName))/255;
        %I = double(imread(caliFileList{3}))/255;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate grid for scatter zoom for QR
        [pY,pSZ] = genIXgrid2(size(I),[800 800],[0 0]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first order snap to center of QR
        pY = fliplr(pY);
        ITER = [2 2];
        [pY] = runZoomSequence(I,2,pY,boxSequence,zoomSequence,map_Q,beta_Q,WINDOW_Q,ITER);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sparkle regress to red corner boxes from end points of first order snap
        fP = pY(:,:,end);
        [sP] = runSparkleSequence(I,2,fP,SPboxSequence,SPzoomSequence,map_S,beta_S,WINDOW_S,[8 2],'','');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % zoom snap to red corners
        ITER = [1 2];
        for centerPoint = 1:size(fP,1)
            for featurePoint = 1:size(sP,1)
                [zoomTrail] = runZoomSequence(I,4,sP(featurePoint,:,centerPoint),boxSequence_F,zoomSequence_F,map_F{featurePoint},beta_F{featurePoint},WINDOW_F{featurePoint},ITER);
                sP(featurePoint,:,centerPoint) = squeeze(zoomTrail(:,:,end));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for good snaps
        msg = {};
        [msg] = readQRdata(I,sP,10,false);
        cnt = 1;
        dy_predict_points = [];
        probeCP = [];
        goodPoints = [];
        GOOD_red_corners = [];
        % find good msg
        for centerPoint = 1:numel(msg)
            if ~isempty(msg{centerPoint})
                goodPoints = [goodPoints;centerPoint];
                probeCP(cnt,:) = fP(centerPoint,:);
                GOOD_red_corners = cat(3,GOOD_red_corners,sP(:,:,centerPoint));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sparkle predict the blue dynamic boxes
                %[dy_predict_points(:,:,cnt)] = runSparkleSequence(I,2,probeCP(cnt,:),DYboxSequence,DYzoomSequence,map_D,beta_D,WINDOW_D,[24*4 2]);
                cnt = cnt + 1;
            end
        end
        %
        qr_oPath = '/mnt/tetra/nate/seedlingDATApile/qrSheets/';
        cropANDsaveQRsheet(I,mean(GOOD_red_corners,3),25,qr_oPath,num2str(e));

        save([qr_oPath num2str(e) '.mat'],'GOOD_red_corners');
    catch
    end
end