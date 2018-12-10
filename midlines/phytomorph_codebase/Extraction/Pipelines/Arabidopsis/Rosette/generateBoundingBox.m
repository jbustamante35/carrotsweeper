function [BB] = generateBoundingBox(I,numBox,autoFlag)
    if ~autoFlag
        for e = 1:numBox
            [J,BB{e}] = imcrop(I);
        end
    end
end

%{
    FilePath = '/home/nate/Downloads/Richard_Pipeline/';
    FileList = {};
    FileExt = {'JPG'};
    FileList = gdig(FilePath,FileList,FileExt,1);

    I = imread(FileList{1});
    BB = generateBoundingBox(I,24,0);
    % fix bounding boxes
    BOX = [];
    for e = 1:numel(BB)
        BOX(e,:) = BB{e}(3:4);
    end
    uB = mean(BOX,1);
    for e = 1:numel(BB)
        BB{e}(3:4) = uB;
    end
    % get subIs
    for e = 1:numel(BB)
        subI(:,:,:,e) = getSquarePots(I,BB{e});
    end
    subI = permute(subI,[1 2 4 3]);
    sz = size(subI);
    d = reshape(subI,[prod(sz(1:3)) sz(4)]);
    subI = ipermute(subI,[1 2 4 3]);
    options = statset('Display','Iter');
    SKIP = 1;
    obj = gmdistribution.fit(double(d(1:SKIP:end,:)),3,'Options',options);
    %% cluster pixels
    for e = 1:size(subI,4)
        cI = classifyImage(double(subI(:,:,:,e)),obj);
        mask = cI == 2;
        mask = logical(bwlarge(mask));
        out = flattenMaskOverlay(subI(:,:,:,e),mask);
        imshow(out,[]);
        waitforbuttonpress
    end
    
    overlayMass(FileList{1},BB,obj,2);
%}