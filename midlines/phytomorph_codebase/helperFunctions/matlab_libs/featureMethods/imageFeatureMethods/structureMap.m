function [lam] = structureMap(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % structureMap --> coded for R2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I       := image
    %           para    := parameters for running the script    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           cM      := corner map       -> corner strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%
    % first order
    D1 = diffmethod(I,para);
    % second degree first order
    D21 = [raster(D1(:,:,1)).^2 raster(D1(:,:,2)).^2 raster(D1(:,:,2)).^2 raster(D1(:,:,1)).^2];
    [lam] = eigenValues(D21);
    lam = reshape(lam,size(D1));
    
    %{
    %%%%%%%%%%%%%%%%%%%%
    % second order
    D2 = diff2method(D1,para);
    D2 = [raster(D2(:,:,1)) raster(D2(:,:,2)) raster(D2(:,:,3)) raster(D2(:,:,4))];
    lam2 = eigenValues(D2);
    lam2 = reshape(lam2,size(D1));
    %}
    
end

%{
    FilePath = '/mnt/spaldingimages/doane/DoaneNSFData/RIL_Data/';
    FileList = {};                                                      % FileList = blanks(1);
    FileExt = {'tif','TIF','bmp','BMP','gif'};    
    verbose = 1;
    fileList = sdig(FilePath,{},FileExt,verbose);




    para.method.value = 'finite';
    para.method.notes = 'differential method';

    
    I = myReader(fileList{1}{1});
    scale = .25;
    [sI BOX] = myCrop(I,scale);

    Ig = rgb2gray(I);
    Ig = double(Ig);
    scale = 1;
    [sIg BOX] = myCrop(Ig,scale,BOX);




scale = 4;
sI = imresize(I,4^-1);
[s BOX] = imcrop(sI,[]);
BOX = BOX * scale;
sIg = imcrop(Ig,BOX);

for s = 5:21
    smooth_sIg = imfilter(double(sIg),fspecial('gaussian',[51 51],s),'replicate');
    lam = structureMap(smooth_sIg,O.op.para.structureMap.para);
    imshow(lam(:,:,2),[]);
    drawnow
end

%}