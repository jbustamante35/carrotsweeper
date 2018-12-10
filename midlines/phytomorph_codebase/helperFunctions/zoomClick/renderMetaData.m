function [msg] = renderMetaData_dev(I,map,beta,WINDOW,map_S,beta_S,WINDOW_S,map_F,beta_F,WINDOW_F,oPath,disp,toDisk)
    try
        mkdir(oPath)
        tic
        if ischar(I)
            [pth,nm,ext] = fileparts(I);
            I = double(imread(I))/255;
        end


        [pY,pSZ] = genIXgrid2(size(I),[800 800],[0 0]);
        pY = fliplr(pY);


        boxSequence = [[1300 1300];[900 900]];
        zoomSequence = [.08 .15];
        ITER = [2 2];
        [pY] = runZoomSequence(I,pY,boxSequence,zoomSequence,map,beta,WINDOW,ITER);

        boxSequence_SP = [1300 1300];
        zoomSequence_SP = [.15];
        fP = pY(:,:,end);
        [sP] = runSparkleSequence(I,fP,boxSequence_SP,zoomSequence_SP,map_S,beta_S,WINDOW_S);


        %
        boxSequence_F = [[300 300];[75 75]];
        zoomSequence_F = [.08 .5];
        ITER = [1 2];
        for centerPoint = 1:size(fP,1)
            for featurePoint = 1:size(sP,1)
                [zoomTrail] = runZoomSequence(I,sP(featurePoint,:,centerPoint),boxSequence_F,zoomSequence_F,map_F{featurePoint},beta_F{featurePoint},WINDOW_F{featurePoint},ITER);
                sP(featurePoint,:,centerPoint) = squeeze(zoomTrail(:,:,end));
            end
        end






        [msg] = readQRdata(I,sP,10,false);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find valid QR reads
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        flag = [];
        for e = 1:numel(msg)
            if ~isempty(msg{e})
                flag(e) = true;
            else
                flag(e) = false;
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find valid QR reads
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = find(flag);
        for e = 1:numel(fidx)


            [qrA,qrW,qrH,whA,whW,whH,huA,huW,huH] = getInformationSheetDimensions(sP(:,:,fidx(e)));


            wholeBOX = [(sP(1,:,fidx(e))) round(whW) round(whH)];
            wholeSheet = mimcrop(I,wholeBOX,round([whW whH]),whA*pi/180);
            %wholeBOX(1:2) = 1;

            dynamicBOX = [(sP(4,:,fidx(e))) round(whW) round(whH)-round(qrH)];
            dynamicSheet = mimcrop(I,dynamicBOX,[round(whW) round(whH)-round(qrH)],whA*pi/180);
            %dynamicBOX(1:2) = (sP(4,:,fidx(e))) - (sP(1,:,fidx(e)));


            qrBOX = [(sP(1,:,fidx(e))) round(qrW) round(qrH)];
            qrSheet = mimcrop(I,qrBOX,round([qrW qrH]),qrA*pi/180);
            %qrBOX(1:2) = 1;

            huBOX = [(sP(2,:,fidx(e))) round(huW) round(huH)];
            huSheet = mimcrop(I,huBOX,round([huW huH]),huA*pi/180);
            huBOX(1:2) = (sP(2,:,fidx(e))) - (sP(1,:,fidx(e)));
            %huBOX(2) = 1;

            %{
            imshow(wholeSheet,[]);
            hold on
            rectangle('Position',qrBOX,'EdgeColor','r');
            rectangle('Position',huBOX,'EdgeColor','g');
            rectangle('Position',dynamicBOX,'EdgeColor','b');
            %}


        end


        if disp
            close all
            imshow(I,[]);
            hold on
            plot(squeeze(pY(:,1,1)),squeeze(pY(:,2,1)),'g*')
            plot(squeeze(pY(:,1,end)),squeeze(pY(:,2,end)),'r*')
            for p = 1:size(pY)
                plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'k')
                plot(squeeze(pY(p,1,:)),squeeze(pY(p,2,:)),'b.')
                hold on
            end
            for p = 1:size(sP,3)
                plot(sP(:,1,p),sP(:,2,p),'k.')
            end
            drawnow
        end
        toc

        for e = 1:numel(msg)
            if ~isempty(msg{e})
                msg{e}
            end
        end



        if toDisk
            imwrite(I,[oPath nm '.tif']);
            save([oPath nm '.mat'],'sP','msg','wholeSheet','fidx')
        end
    catch ME
        ME
    end
end