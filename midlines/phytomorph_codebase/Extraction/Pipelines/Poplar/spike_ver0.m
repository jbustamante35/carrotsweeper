fileName = '/home/nate/Downloads/SUCCESS #7 G-GECO root + GSEs 20+ min timecourse.czi'; 
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
    %%
    fs = fs / (2^12-1)
    %%
    close all
    for e = 1:size(fs,4)
        imshow(fs(:,:,:,e),[]);
        drawnow
    end
    %%
    sMap= std(squeeze(fs(:,:,2,:)),1,3);
    %%
    sMap= std(squeeze(fs(:,:,2,:).*fs(:,:,1,:)),1,3);
    M = sMap > graythresh(sMap)*.7;
     M = bwareaopen(M,100);
     %M = imdilate(M,strel('disk',15,0));
    %%
    imshow(M,[]);
    %%
    close all
    for e = 1:size(fs,4)
        out = flattenMaskOverlay(fs(:,:,:,e),M);
        imshow(out,[]);
        drawnow
    end
    %%
    writerObj = VideoWriter('spike1.avi');
    open(writerObj);
    close all
    R = regionprops(M,'PixelIdxList');
    CL = {'r' 'b' 'r' 'c' 'm' 'w' 'y'};
    for e = 1:size(fs,4)
        out = fs(:,:,:,e);
        E = {};
        for r = 1:numel(R)
            Z = zeros(size(fs,1),size(fs,2));
            Z(R(r).PixelIdxList) = 1;
            Z = imdilate(Z,strel('disk',21,0));
            E{r} = bwboundaries(Z);
            idx = find(Z);
            tmp = fs(:,:,1,e);
            V(e,r) = mean(tmp(idx));
            
            
            
            
            out = flattenMaskOverlay(out,logical(Z),.25,CL{r});
            
        end
        %plot(V)
        %drawnow
        
        imshow(out,[]);
        hold on
        for r = 1:numel(E)
            plot(E{r}{1}(:,2),E{r}{1}(:,1),'b');
        end
        hold off
            drawnow
        frame = getframe;
        writeVideo(writerObj,frame);
        
        
    end
    close(writerObj);