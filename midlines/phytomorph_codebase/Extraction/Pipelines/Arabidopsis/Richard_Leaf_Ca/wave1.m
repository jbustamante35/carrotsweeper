function [] = wave1(fileName,oPath,disp)
    %% read a file
    fprintf(['Start loading data.\n']);
    %fileName = '/mnt/snapper/nate/forRichard/forACA/160529_aca4.1-11.1 T1 plantE.nd2';
    meta = imreadBFmeta(fileName);
    [pth,nm,ext] = fileparts(fileName);
    MF = meta.nframes;
    %MF = 50;
    D = imreadBF(fileName,1,[1:MF],1);
    E = imreadBF(fileName,1,[1:MF],2);
    %replace zeros with 1
    D(D(:)==0) = 1;
    % make ratio
    R = E.*D.^-1; 
    fprintf(['End loading data.\n']);
    if disp
        %% view ratio movie
        close all
        for e = 1:size(R,3)
            imshow(R(:,:,e),[1.2 3.2]);
            drawnow
        end
    end
    %% smooth in time
    fprintf(['Start time domain analysis.\n']);
    % smooth each pixel value over time
    fprintf(['Smooth data sig over time.\n']);
    tsR = timeSmooth(R,11);
    % take the 1st derivative of the time signal
    fprintf(['Take 1st der.\n']);
    dtsR = DGradient(tsR,1,3);
    % smooth the 1st derivative
    fprintf(['Smooth 1st der.\n']);
    sdtsR = timeSmooth(dtsR,11);
    % take the 1st derivative of the 1st derivative = 2nd
    fprintf(['Take 2nd der.\n']);
    ddR = DGradient(sdtsR,1,3);
    % smooth the 2nd derivative
    fprintf(['Smooth 2nd der.\n']);
    sddR = timeSmooth(ddR,11);
    fprintf(['End time domain analysis.\n']);
    %% smoothing dJ
    fprintf(['Start space domain analysis.\n']);
    ssdJ = zeros(size(sdtsR));
    for e = 1:size(dtsR,3)
        fprintf(['Start the:' num2str(e) 'frame.\n']);
        ssdJ(:,:,e) = imfilter(sdtsR(:,:,e),fspecial('disk',21),'replicate');
        fprintf(['End the:' num2str(e) 'frame.\n']);
    end
    % smoothing d2J
    ssd2J = zeros(size(sddR));
    for e = 1:size(dtsR,3)
        fprintf(['Start the:' num2str(e) 'frame.\n']);
        ssd2J(:,:,e) = imfilter(sddR(:,:,e),fspecial('disk',21),'replicate');
        fprintf(['End the:' num2str(e) 'frame.\n']);
    end
    fprintf(['End space domain analysis.\n']);
    %% find peaks for first der
    fprintf(['Start find peaks of first derivative.\n']);
    sig = reshape(tsR,[size(tsR,1)*size(tsR,2) size(tsR,3)]);
    dsig = imdilate(sig,ones(1,61));
    front0 = dsig == sig;
    front0 = reshape(front0,size(sdtsR));
    fprintf(['End find peaks of first derivative.\n']);
    %% find peaks for first der
    fprintf(['Start find peaks of first derivative.\n']);
    sig = reshape(ssdJ,[size(sdtsR,1)*size(sdtsR,2) size(sdtsR,3)]);
    dsig = imdilate(sig,ones(1,61));
    front = dsig == sig;
    front = reshape(front,size(sdtsR));
    fprintf(['End find peaks of first derivative.\n']);
    %% find peaks for second der
    fprintf(['Start find peaks of second derivative.\n']);
    sig = reshape(ssd2J,[size(sdtsR,1)*size(sdtsR,2) size(sdtsR,3)]);
    dsig = imdilate(sig,ones(1,61));
    front2 = dsig == sig;
    front2 = reshape(front2,size(sdtsR));
    fprintf(['End find peaks of second derivative.\n']);
    %% get connected events
    BLOB = ssdJ > .006;
    CC = bwconncomp(BLOB);
    %% filter front into wave front
    wfront = front.*BLOB;
    CCf = bwconncomp(wfront);
    %% filter front2 
    BLOB2 = ssd2J > .5*10^-3;
    wfront2 = front2.*BLOB2;
    %% filter front2 
    BLOB0 = tsR > 2;
    wfront0 = front0.*BLOB0;
    %%

    %% show front1 and front2
    close all
    %MOV = zeros(size(tsR,1),size(tsR,2),3,size(tsR,4));
    v = VideoWriter([oPath nm '.avi']);
    open(v)
    for e = 1:size(front2,3)
        tmp = tsR(:,:,e);
        tmp = imfilter(tmp,fspecial('disk',11));
        out = flattenMaskOverlay(bindVec(tmp),logical(wfront(:,:,e)),.5,'b');
        out = flattenMaskOverlay(out,logical(wfront2(:,:,e)),.5,'r');
        out = flattenMaskOverlay(out,logical(wfront0(:,:,e)),.5,'g');
        imshow(out,[]);
        title(num2str(e))
        drawnow
        writeVideo(v,out);
        %MOV(:,:,:,e) = out;
    end
    close(v);
end