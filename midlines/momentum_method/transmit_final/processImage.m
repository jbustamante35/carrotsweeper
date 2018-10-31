function [out] = processImage(fileName,disp)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove the background
        I = double(imread(fileName));
        if strcmp(class(I),'uint8')
            I = I / 255;
        elseif strcmp(class(I),'uint16')
            I = I / (2^16-1);
            I = imcomplement(I);
        end
        
        if size(I,3) == 3
            I = rgb2gray(I);
        end
        I(:,end) = [];
        I = handleFLIP(I,[]);
        % added to remove the bad pixel columns in image.
        I(:,1) = [];
        
        Io = I;
        Is = imresize(I,.25);
        BK = imclose(Is,strel('disk',21));
        BK = imfilter(BK,fspecial('disk',21),'replicate');
        BK = imresize(BK,size(Io));
        I = I - BK;
        I = I - min(I(:));
        I = I / max(I(:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain the binary image
        I = imfilter(I,fspecial('gaussian',21,6),'replicate');        
        thresh = graythresh(I);
        skel = I < thresh;
        skel = bwareaopen(skel,500);
        cB = imclearborder(skel);
        skel = skel & ~ cB;
        [g1 g2] = gradient(I);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the contours
        curve = getLevelContours(double(skel),[1 1]);
        ridx = find([curve.length] < 100);
        curve(ridx) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % find the highest curvature for each contour
        KSNIP = 50;
        SMOOTH_VALUE = 20;
        for e = 1:numel(curve)
            o = cwtK(curve(e).data',{SMOOTH_VALUE});
            [J tipIDX{e}] = min(o.K);
            o = cwtK(curve(e).data',{20});
            [J fine_tipIDX] = min((o.K(tipIDX{e}-KSNIP:tipIDX{e}+KSNIP)));
            tipIDX{e} = tipIDX{e} + (fine_tipIDX - KSNIP - 1);
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make distance transform
        B = double(bwdist(~skel));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trace each midline
        for e = 1:numel(curve)
            d1X1 = cwt(curve(e).data(1,:),5,'gaus1');
            d1X2 = cwt(curve(e).data(2,:),5,'gaus1');
            t1 = -d1X2(tipIDX{e});
            t2 = d1X1(tipIDX{e});
            T = [t1 t2];
            T = T / norm(T);
            N = [T(2) -T(1)];
            initD = [T;N];                             %(I,P,initD,maxStep,RHO,RAD,pointDensity,wsigma)
            midlines(e).data = trackFromPointAtGradient(B,curve(e).data(:,tipIDX{e}),initD,10000,20,pi/2,[20 200],.7);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display
        if ~isempty(disp)
            figure(disp)
            hold off
            imshow(Io,[]);
            hold on
            for e = 1:numel(curve)
                plot(curve(e).data(1,:),curve(e).data(2,:),'r');
                plot(midlines(e).data(1,:),midlines(e).data(2,:),'b');
            end

            for e = 1:numel(curve)
                plot(curve(e).data(1,tipIDX{e}),curve(e).data(2,tipIDX{e}),'g*');
            end
            drawnow
        end
        out.midlines = midlines;
        out.contours = curve;
        
        
    catch ME;
        out.ME = ME;
    end
end