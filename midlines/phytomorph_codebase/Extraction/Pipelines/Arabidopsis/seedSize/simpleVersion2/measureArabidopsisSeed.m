function [] = measureArabidopsisSeed(fileName,oPath)
    % make directory
    fprintf(['Make output directory.\n']);
    mkdir(oPath);
    fprintf(['Read Image.\n']);
    if ischar(fileName)
        I = imread(fileName);
    else
        I = fileName;
    end
    fprintf(['Transform RGB-->HSV.\n']);
    H = rgb2hsv(I);
    fprintf(['Threshold image.\n']);
    fprintf(['Image size:' num2str(size(H)) '\n']);
    M = H(:,:,3) > graythresh(H(:,:,3));
    M = bwareaopen(M,50);
    fprintf(['Measure with region props2.\n']);
    R = regionprops(M,'Area','PixelIdxList','MajorAxisLength','MinorAxisLength');
    cidx = count([R.Area]);
    cidx1 = count([R.MajorAxisLength]);
    cidx2 = count([R.MinorAxisLength]);
    fidx = find(cidx==1 & cidx1==1 & cidx2==1);
    singleMask = zeros(size(M));
    for e = 1:numel(fidx)
        singleMask(R(fidx(e)).PixelIdxList) = 1;
    end
    fprintf(['Make overlay.\n']);
    out = flattenMaskOverlay(I,logical(singleMask),.6,'r');
    fprintf(['Get file parts.\n']);
    if ischar(fileName)
        [p,nm,ext] = fileparts(fileName);
    else
        nm = 'test';
    end
    fprintf(['Write image data.\n']);
    imwrite(out,[oPath nm '_return.jpg']);
    phenotypeData = [[R(fidx).Area]' [R(fidx).MajorAxisLength]' [R(fidx).MinorAxisLength]'];
    fprintf(['Write Data.\n']);
    csvwrite([oPath nm '_phenotypeData.csv'],phenotypeData);
end

%{
    fileName = '/home/nate/10_26_trials/GS_test_COL_4800dpi_rgb028.tif';
    oPath = '/home/nate/10_26_trials/';
    measureArabidopsisSeed(fileName,oPath);

    I = imread('/home/nate/10_26_trials/GS_test_COL_4800dpi_rgb028_return.jpg');

    FilePath = '/home/nate/10_26_trials/Calibration Images/';
    oPath = '/home/nate/10_26_trials/Calibration Images_return/';
    mkdir(oPath);
    FileList = {};
    FileExt = {'tiff','TIF'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    for e = 1:numel(FileList)
        measureArabidopsisSeed(FileList{e},oPath);
    end

    fileName = '/iplant/home/turnersarahd/Clement_drought_scans/Normal/035255243056';
    measureArabidopsisSeed(fileName,'');

%}


