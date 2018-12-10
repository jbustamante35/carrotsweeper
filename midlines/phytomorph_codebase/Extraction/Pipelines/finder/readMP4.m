function [mov] = readMP4(fileName,skip)
    videoR = VideoReader(fileName);
    tI = read(videoR,1);
    % nFrames = 800; first movie
    nFrames = 680; % second movie
    videoR = VideoReader(fileName);
   
    cVec = 1:skip:(nFrames);
    mov = zeros([size(tI) numel(cVec)]);
    cnt = 1;
    for n = 1:numel(cVec)
        mov(:,:,:,cnt) = read(videoR,cVec(n));
        cnt = cnt + 1
    end
end