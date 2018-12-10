function [T ufT BB PS MT sM] = measureKernelLength(I,numberCobs,RAD,gridSites,defaultAreaPix,fill,CHUNK)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                measureKernelLength.m is main function to handle cob analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box and color circle. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                rgb2hsv_fast.m, measureImage.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                I:              The blue header remvoed image to be analyzed.
                numberCobs:     Number of cobs that are expected to be analyzed. 
                RAD:            The value for window size.
                gridSites:      The number of down sample grid sites.
                defaultAreaPix: The default pixel to be considered noise relative to 1200 dpi.
                fill:           The radius of disk for Kernel of an image close operation.
                CHUNK:          The number of chunk for input for FFT in myBlock0.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    %%%%%%%%%%%%%%%%%%%%%%%
    % init return vars 
    ufT = [];
    PS = [];
    PS.widthProfile = [];
    sM = {};
    %%%%%%%%%%%%%%%%%%%%%%%
    try        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % get hsv slice, filter and mask
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting mask creation. \n']);
        %h = fspecial('gaussian',[31 31],11);
        %I = imfilter(I,h);
        fI = rgb2hsv_fast(I,'single','V');
        level = graythresh(fI);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % get binary, close, fill holes, remove small objects
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        B = fI > level;
        B = imopen(B,strel('disk',fill));
        B = imclose(B,strel('disk',fill));
        B = imfill(B,'holes');
        % changed to work with QR and with smaller DPI Dec. 28 2017
        %B = bwareaopen(B,defaultAreaPix);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Feb 27,2018
        % switched to clear border first-->find large
        % find the N largest ears
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        B = imclearborder(B);
        B = bwlarge(B,numberCobs);
        R = regionprops(B,'PixelIdxList','PixelList','Area','Image','BoundingBox');
        fprintf(['Found ' num2str(numel(R)) ' ears \n']);
        %{
        % added for QR codes - start
        if numel(R) > numberCobs
            [~,sidx] = sort([R.Area]);
            R(sidx(1:(end-numberCobs))) = [];
        end
        numberCobs = min(numberCobs,numel(R));
        %}
        % added for QR codes - end
        fprintf(['ending mask creation. \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % iterate over the cobs        
        for e = 1:numberCobs
            fprintf(['starting with ear ' num2str(e) ':' num2str(numberCobs) '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get bounding box - color
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BB{e} = R(e).BoundingBox;
            subI = R(e).Image;
            grayImage = single(imcrop(I,R(e).BoundingBox))/255;
            grayImage = rgb2gray(grayImage);
            fprintf(['Operating on sub image of size: [' num2str(size(grayImage)) ']:[' num2str(size(I)) ']\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get width profile
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            WIDTH = sum(subI,2);
            fidx = find(WIDTH ~= 0);
            uW = mean(WIDTH(fidx));
            PS.average_WIDTH(e) = uW;
            % store width profile
            PS.widthProfile = [PS.widthProfile ;interp1(1:numel(WIDTH),WIDTH,linspace(1,numel(WIDTH),1000))];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display
            % process the grayscale image, take gradient, look at pos and neg
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = fspecial('average',41);            
            grayImage = imfilter(grayImage,h);
            % look at gradient of grayscale image - could be yellow
            [~,grayImage] = gradient(grayImage);
            % smooth cols gradient
            grayImage = imfilter(single(grayImage),fspecial('gaussian',[41 41],11));
            % find the beginning and ending of kernels - pos and neg
            gMSK = grayImage > 0;
            lMSK = grayImage < 0;
            % strip the image border to zero
            subI(1,:) = 0;
            subI(:,end) = 0;
            subI(end,:) = 0;
            % init vars to measuure image
            tG = {};
            tL = {};
            parfor r = 1:numel(RAD)
                fprintf(['starting with fft window ' num2str(r) ':' num2str(numel(RAD)) '\n']);
                % set the current window size
                dR = RAD(r);
                % errode such that the fft window samples only ear image
                toMeasure = imerode(subI,strel('rectangle',[2*dR+1 20]));
                % call to measure kernel period of the rising edge
                [Tg tG{r}] = measureImage(gMSK.*grayImage,toMeasure,gridSites,dR,CHUNK);
                % call to measure the kernel period of the falling edge
                [Tl tL{r}] = measureImage(lMSK.*grayImage,toMeasure,gridSites,dR,CHUNK);
                % stack results together
                MT(r,:,e) = [Tg Tl];
                % average of rising and falling edge period
                T(e,r) = nanmean([Tg Tl]);                
                fprintf(['ending with fft window ' num2str(r) ':' num2str(numel(RAD)) '\n']);
            end
            fprintf(['ending with ear ' num2str(e) ':' num2str(numberCobs) '\n']);
            ufT.G{e} = tG;
            ufT.L{e} = tL;
        end
    catch ME
        getReport(ME)
        fprintf(['******error in:measureKernelLength.m******\n']);
    end
end