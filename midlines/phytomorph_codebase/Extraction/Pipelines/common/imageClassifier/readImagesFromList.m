function [imageStack] = readImagesFromList(FileList,reducPer)
    numToRead = numel(FileList);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read training image stack - X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['************************************************************************\n']); 
    fprintf(['************************************************************************\n']);
    fprintf(['************************************\n']);
    fprintf(['start reading images \n']);
    fprintf(['************************************\n']);
    tmpI = imread(FileList{1});
    tmpI = imresize(tmpI,reducPer);
    sz = size(tmpI);
    if numel(sz) == 2
        sz = [sz 1];
    end
    imageStack = zeros([sz numToRead]);
    for e = 1:numToRead
        fprintf(['Start reading image:' num2str(e) '\n']);
        tmpI = imread(FileList{e});
        tmpI = imresize(tmpI,reducPer);
        imageStack(:,:,:,e) = tmpI;
        fprintf(['Stop reading image:' num2str(e) '\n']);
    end
    fprintf(['************************************\n']);
    fprintf(['stop reading images \n']);
    fprintf(['************************************\n']);
    fprintf(['************************************************************************\n']); 
    fprintf(['************************************************************************\n']);
end