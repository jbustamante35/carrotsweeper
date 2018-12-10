function [] = cellForDirk(fileName,oPath)
    mkdir(oPath);
    info = imfinfo(fileName)
    for e = 1:numel(info)
        I(:,:,:,e) = imread(fileName,e);
    end
    U = double(mean(I,4))/(2^12-1);
    M = double(max(I,[],4))/(2^12-1);
    F = mean(cat(4,.05*M,.95*U),4);
    
    E = edge(F,'Canny',.05);
    E = imclose(E,strel('disk',2,0));
    E = bwareaopen(E,100,8); 
    C = ~E;
    innerM = bwareaopen(C,200,4);
    close all
    innerM2 = (innerM==0 & C == 1);
    T = cat(3,innerM,C,innerM2);
    C = C | innerM;
    L = bwlabel(C,4);
    RGB = label2rgb(L);
    R = regionprops(L,'PixelIdxList','Area');
    Area = [R.Area];
    csvwrite([oPath 'area.csv'],Area');
    imwrite(RGB,[oPath 'image.tif']);
end