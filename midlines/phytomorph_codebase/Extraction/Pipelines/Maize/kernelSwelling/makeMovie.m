inFilePath = '/mnt/snapper/kernelSwellingData/Scott/rawData/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);

%% sort
for e = 1:numel(SET)
    NM = [];
    for img = 1:numel(SET{e})
        [p,nm,ext] = fileparts(SET{e}{img});
        NM(img) = str2num(nm);
    end
    [~,sidx] = sort(NM);
    SET{e} = SET{e}(sidx);
end
%%
I = imread(SET{1}{1});
[J,BOX] = imcrop(I);

%% stack frames
close all
SKIP = 10;
cnt = 1;
P = zeros(2*348,3*339,3,390);
writerObj = VideoWriter('swell_mass2.avi');
writerObj.FrameRate = 4;
open(writerObj);
for img = 1:SKIP:390
    tmpP1 = [];
    tmpP2 = [];
    tmpP3 = [];
    tmp = {};
    OFFSET = 6;
    parfor e = 1:6
        tmp{e} = imread(SET{e+OFFSET}{img},'PixelRegion',{[1 30 10440], [1 30 10164]});        
    end
    
    
    for e = 1:6    
        if e == 1 || e == 2
            tmpP1 = cat(1,tmpP1,tmp{e});
        elseif e == 5 || e == 6
            tmpP2 = cat(1,tmpP2,tmp{e});
        elseif e == 3 || e == 4
            tmpP3 = cat(1,tmpP3,tmp{e});
        end
        e
    end
    
    
    
    P(:,:,:,cnt) = cat(2,tmpP1,tmpP2,tmpP3);    
    imshow(P(:,:,:,cnt)/255,[],'Border','tight');    
    drawnow
    cnt = cnt + 1;
    
    frame = getframe;
    writeVideo(writerObj,frame);
    img
end
close(writerObj);
%%  make movie again
writerObj = VideoWriter('swell_mass.avi');
writerObj.FrameRate = 3;
open(writerObj);
for t = 1:size(P,4)
    imshow(P(:,:,:,t)/255,[],'Border','tight')
    drawnow    
    frame = getframe;
    writeVideo(writerObj,frame);    
end
close(writerObj);
%%
mov = immovie(P);
%%
movie2avi(mov, 'myPeaks.avi', 'compression', 'None')