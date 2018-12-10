function [toP K] = confocalRatio(fileName,toDISK,oPath,disp)
    if isdeployed
        javaaddpath([pwd filesep 'bioformats_package.jar']);
    end

    % read the image stack
    fs = [];
    data = bfopen(fileName);
    [pth,nm,ext] = fileparts(fileName);
    for e = 1:(size(data{1},1)/3)
        idxB = (e-1)*3;
        for k = 1:3
            fs(:,:,k,e) = fliplr(data{1}{idxB+k,1});
        end
    end
    fs = double(fs);
    sig = mean(mean(fs,1),4);
    sig = squeeze(sig);
    sig = sig(:,1);
    sig = sig / sum(sig);
    PAD = 100;
    [J,sidx] = max([mean(sig(1:PAD)),mean(sig(end-PAD:end))]);
    if sidx == 2
        fs = flipdim(fs,2);
    end
    if disp
        h1 = figure;
        h2 = figure;
    end
    % process image stack
    toP = [];
    K = [];
    tMASK = [];
    for e = 1:size(fs,4)-5
        tm = clock;
        % pop image from stack
        sourceO = fs(:,:,:,e)/255;
        % mean over channels 1:2
        source = mean(sourceO(:,:,1:2),3);
        % make mask
        MASK = source > .8*graythresh(source);
        % area open
        MASK = bwareaopen(MASK,1000);
        % sum mask
        tipMask = sum(MASK,1);
        % tip mask stack
        tMASK = [tMASK,tipMask'];
        % average vertical signal
        sig1 = mean(source,1);
        % filter vertical signal
        sig1 = imfilter(sig1,fspecial('disk',11),'replicate');
        % second filter
        sig2 = imfilter(sig1,fspecial('disk',51),'replicate');
        % gradient for super smooth signal
        gsig = gradient(sig2);
        % find max of super smooth signal
        [~,gidx] = max(sig2);
        % find max of gradient signal
        [~,gidx2] = max(gsig);
        % obtain threshold
        uT = mean(sig1((end-30):end));
        % magnify threshold
        uT = uT*1.2;
        % threshold signal
        th1 = sig1 < uT;
        th1 = tipMask < 5;
        % fill from base of root to max signal location
        th1(1:gidx) = 0;
        % create heaveside functions for base to max peak
        th2 = zeros(size(th1));
        % create heaveside function for base to max grad
        th3 = zeros(size(th1));
        % fill heaveside
        th2(1:gidx) = 1;
        % fill HS
        th3(1:gidx2) = 1;
        % find the locations for whole rooot
        fidx = find(th1);
        % reset whole root mask
        th1 = zeros(size(th1));
        % fill root mask 
        th1(1:fidx(1)) = 1;
        % make square mask
        squareMASK = zeros(size(source));
        % fillsquare mask
        squareMASK(:,gidx2:fidx(1)) = 1;
        % multiply the square mask by the root mask
        MASK = MASK.*squareMASK;
        % filter the soure image
        sourceO = imfilter(sourceO,fspecial('disk',5),'replicate');
        % make the ratio image
        RATIO = sourceO(:,:,2).*sourceO(:,:,1).^-1;
        % mask the ratio image
        toM = MASK.*RATIO;
        % create mesh strip
        Y = sum(toM,1);
        X = sum(MASK,1);
        Y = Y.*X.^-1;
        % clip mesh strip and interp the mesh strip
        xidx = find(X);
        Y = Y(xidx);
        k = interp1(Y,linspace(1,numel(Y),200));
        % store the mesh strip
        K(:,e) = k';
        % average over the tip
        toP(e) = sum(toM(:))/sum(MASK(:));
        
        if disp
            mag = 100;
            % create overlay
            out = flattenMaskOverlay(source,logical(MASK));
            figure(h1);
            imshow(out,[])
            hold on
            plot(sig2*mag*10,'y')
            plot(th1*mag/10,'r')
            plot(th2*mag/10,'b');
            plot(th3*mag/10,'c');

            hold off
            figure(h2);
            plot(toP);
            drawnow
        end
        
    end
    l = [];
    for th = 5:30
        for e = 1:size(tMASK,2)
            fidx = find(tMASK(:,e) > th);
            l(th-4,e) = fidx(end);
        end
    end
    dL = diff(l,1,2);
    [~,midx] = min(sum(dL.*dL,2));
    if toDISK
        mkdir(oPath);
        csvwrite([oPath nm '_rootTipSignal.csv'],toP);
        csvwrite([oPath nm '_kymogrpah.csv'],K);
        csvwrite([oPath nm '_tipProfile.csv'],tipMask);
        csvwrite([oPath nm '_growthRate.csv'],l(midx,:));
    end
end