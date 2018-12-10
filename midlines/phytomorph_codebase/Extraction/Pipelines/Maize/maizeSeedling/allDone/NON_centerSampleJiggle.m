function [wXi,wYi] = NON_centerSampleJiggle(trainTable,DELTA,scaleRange,OVER_SQUARE,reSZ2,totalR)

    % sample stores
    wX = {};
    wY = {};
    parfor e = 1:size(trainTable,1)
        tic;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image
        I = double(imread(trainTable.wholeImageFileName{e}))/255;
        
        
        [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,20);
        
        MSK = zeros(size(I,1),size(I,2));
        MSK2 = zeros(size(I,1),size(I,2));
        MSK(1:round(cropLine),:) = 1;
        MSK(:,round(qrCropBox(1)):round(qrCropBox(1)+qrCropBox(3))) = 0;
        MSK2(round(cropLine):end,:) = 1;
        
        
        % snip of side of purple
        I(:,(end-50):end,:) = [];
        % get meta data box
        qrBOX = trainTable.metaDataBox{e};
        % center point of the box
        CP_BOX = qrBOX(1:2) + .5*qrBOX(3:4);
        % center point of the image
        CP = size(I)/2;
        % remove the color channel size
        CP(3) = [];
        % flip dim forwork
        CP = flipdim(CP,2);
        % shift to center parameter
        SHIFT = round(CP - CP_BOX);
        % init the stores for r samples
        LOCY = [];
        subI = [];
        % shift the image to the center
        shI = circshift(I,SHIFT(1),2);
        shI = circshift(shI,SHIFT(2),1);
        
        % shift the image to the center
        MSK = circshift(MSK,SHIFT(1),2);
        MSK = circshift(MSK,SHIFT(2),1);
        
        MSK2 = circshift(MSK2,SHIFT(1),2);
        MSK2 = circshift(MSK2,SHIFT(2),1);
        
        %fidx = find(MSK);
        %fidx2 = find(MSK2);
        fidx = [];
        fidx2 = [];
        [fidx(:,2),fidx(:,1)] = find(MSK);
        [fidx2(:,2),fidx2(:,1)] = find(MSK2);
        
        %{
        MSK = zeros(size(shI,1),size(shI,2));
        MSK(CP(2)-600:CP(2)+600,CP(1)-600:CP(1)+600) = 1;
        DIST = bwdist(MSK);
        DIST(:,1:CP(1)) = 0;
        fidx = [];
        [fidx(:,2),fidx(:,1)] = find(DIST > 30);
        %}
        
        
        
        
        % generate some number of repeats per image
        cnt = 1;
        for r = 1:totalR
            % generate a random scale
            meanScale = mean(scaleRange);
            rng = scaleRange(2) - scaleRange(1);
            tmpScale = (rng*(rand(1) - .5) + meanScale);
            ridx = randi(size(fidx,1),1);
            ridx2 = randi(size(fidx2,1),1);
       
            CPs = fidx(ridx,:);
            CPs2 = fidx2(ridx2,:);
            
            
            % scale the crop box
            oBOX = point2Box((CPs),round(OVER_SQUARE*tmpScale^-1));
            % crop the image
            tmp = mimcrop(shI,oBOX,reSZ2);
            
            % scale the crop box
            oBOX2 = point2Box((CPs2),round(OVER_SQUARE*tmpScale^-1));
            % crop the image
            tmp2 = mimcrop(shI,oBOX2,reSZ2);
            
            
            
            % store the image
            subI(:,:,:,cnt) = tmp;
            cnt = cnt + 1;
            
            %{
            % for view
            imshow(tmp,[]);
            drawnow
            pause(.5)
            %}
            
            subI(:,:,:,cnt) = tmp2;
            cnt = cnt + 1;
            
            %{
            % for view
            imshow(tmp2,[]);
            drawnow
            pause(.5)
            %}
            
        end
        e
        wX{e} = subI;
        wY{e} = LOCY;
        fprintf(['Total time is: ' num2str(toc*size(trainTable,1)/60/12) '\n']);
    end

    % stack the QR center location data
    sz = size(wX{1});
    per = sz(4);
    sz(4) = sz(4)*numel(wX);
    wXi = zeros(sz);
    wYi = [];
    for e = 1:numel(wX)
        str = (e-1)*per + 1;
        stp = str + per - 1;
        wXi(:,:,:,str:stp) = wX{e};
        tmp = wY{e};
        sz = size(tmp);
        %sz(3) = 1;
        %tmp = reshape(tmp,[(sz(1:2)) sz(3)]);
        wYi = [wYi;tmp];
    end
end
