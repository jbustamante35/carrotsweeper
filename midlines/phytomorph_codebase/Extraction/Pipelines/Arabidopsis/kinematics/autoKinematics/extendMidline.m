function [midlineM,rootWidth] = extendMidline(I,path)
    % smooth the midline
    path = imfilter(path,fspecial('average',[1 100]),'replicate');
    path = arcLength(path','arcLen')';

    
    
    
    WIDTH_NUMP = 600;
    PCA_RHO = 50;
    WIDTH = 300;
    Domain = genCurvilinearDomain(path',PCA_RHO,WIDTH,WIDTH_NUMP,[]);

    dX = -diff(Domain,1,1);
    SNIP = 200;
    dX = mean(dX(1:SNIP,:,:),1);
    dNOR = sum(dX.^2,3).^-.5;
    dX = bsxfun(@times,dX,dNOR);
    EXT = 500;
    EXT = linspace(0,EXT,EXT);
    EXT = bsxfun(@times,EXT',dX);
    EXT = bsxfun(@plus,EXT,Domain(1,:,:));
    
    
    mDomain = cat(1,flipdim(EXT,1),Domain);
    disp = 0;
    if disp
        SKIP = 10;
        imshow(I,[]);
        hold on
        for p = 1:SKIP:size(Domain,1)
            plot(Domain(p,:,1),Domain(p,:,2),'r')
        end
        for p = 1:SKIP:size(EXT,1)
            plot(EXT(p,:,1),EXT(p,:,2),'b')
        end
    end
    
    % make mask for kinematics
    F1 = mDomain(:,:,1);
    F2 = mDomain(:,:,2);
    FT = ba_interp2(double(I)/255,F1,F2);
    WIDTH = mean(FT,1);
    LENGTH = mean(FT,2);
    sL = std(FT,1,2);
    thresh = graythresh(WIDTH);
    thresh2 = graythresh(sL);
    MSK = double(sL > thresh2) * double(WIDTH < thresh);
    MSK = cumsum(MSK,1) > 1;
    R = regionprops(logical(MSK),'Area','PixelIdxList');
    [J,midx] = max([R.Area]);
    MSK = zeros(size(MSK));
    MSK(R(midx).PixelIdxList) = 1;
    
    
    
    
    [midlineM{1}] = midRib_ver0(FT,MSK,F1,F2);
    [midlineM{2}] = midRib_ver1(FT,MSK,F1,F2);
    [midlineM{3}] = midRib_ver3(FT,MSK,F1,F2);
    WID = find(any(MSK,1));
    rootWidth = mean(WID);
    
    
    
    % trace box mask
    dB = bwboundaries(MSK);
    dB = dB{1};
    idx = sub2ind([size(mDomain,1) size(mDomain,2)],dB(:,1),dB(:,2));
    BOX = [];
    BOX(:,1) = F1(idx);
    BOX(:,2) = F2(idx);

end