function [D P] = clusterMovie(videoFile,GMM,reSize,disp)
    readerObj = VideoReader(videoFile);
    cnt = 1;
    while readerObj.hasFrame
        tmp = double(readerObj.readFrame)/255;
        %tmp = imfilter(tmp,fspecial('gaussian',[21 21],7),'replicate');
        if reSize ~= 1
            tmp = imresize(tmp,reSize);
        end
        sz = size(tmp);
        tmp = reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]);
        [tmp NL] = cluster(GMM,tmp);
        tmp = reshape(tmp,sz(1:2));
        if reSize ~= 1
            tmp = imresize(tmp,reSize.^-1);
        end
        
        %D{cnt} = sparse(tmp==4);
        
        D(:,:,:,cnt) = tmp;
        P(:,:,:,cnt) = NL;
        cnt = cnt +1;
        if disp
            tmp = label2rgb(tmp);
            imshow(tmp,[])
            drawnow
        end
        cnt
    end
end