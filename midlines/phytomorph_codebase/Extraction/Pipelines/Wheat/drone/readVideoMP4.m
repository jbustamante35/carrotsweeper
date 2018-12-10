function [D] = readVideoMP4(videoFile,reSize)
    readerObj = VideoReader(videoFile);
    cnt = 1;
    while readerObj.hasFrame
        tmp = double(readerObj.readFrame)/255;
        if reSize ~= 1
            tmp = imresize(tmp,reSize);
        end
        D(:,:,:,cnt) = tmp;
        cnt = cnt +1;
    end
end