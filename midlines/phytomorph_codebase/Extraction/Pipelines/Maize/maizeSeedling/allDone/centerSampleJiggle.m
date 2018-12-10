function [wXi,wYi] = centerSampleJiggle(trainTable,DELTA,scaleRange,OVER_SQUARE,reSZ2,totalR)

    % sample stores
    wX = {};
    wY = {};
    parfor e = 1:size(trainTable,1)
        tic;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image
        I = double(imread(trainTable.wholeImageFileName{e}))/255;
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
        % generate some number of repeats per image
        for r = 1:totalR
            % generate a random scale
            meanScale = mean(scaleRange);
            rng = scaleRange(2) - scaleRange(1);
            tmpScale = (rng*(rand(1) - .5) + meanScale);
            fprintf(['temp scale is:' num2str(tmpScale) '\n'])
            % generate random jiggle
            JG = (DELTA*(rand([1 2]) - .5));
            JG = round(JG);
            % shift the image to the center
            shI = circshift(I,SHIFT(1)+JG(1),2);
            shI = circshift(shI,SHIFT(2)+JG(2),1);
            % scale the crop box
            oBOX = point2Box((CP),round(OVER_SQUARE*tmpScale^-1));
            % store the values
            LOCY(r,:) = [JG*tmpScale tmpScale];
            % crop the image
            tmp = mimcrop(shI,oBOX,reSZ2);
            % store the image
            subI(:,:,:,r) = tmp;
            % for view
            %imshow(tmp,[]);
            %drawnow
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
