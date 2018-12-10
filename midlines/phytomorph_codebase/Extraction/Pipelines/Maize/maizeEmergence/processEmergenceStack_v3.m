function [miniStack,miniMask] = processEmergenceStack_v3(imageStack,oPath)
    fprintf(['****************************************************************************************\n']);
    versionString = ['Starting emergence analysis algorithm. \nPublication Version 1.0 - Monday, April 3, 2017. \n'];
    fprintf(versionString);
    fprintf(['****************************************************************************************\n']);
    
    % constant to tell the number of frames to average over when finding
    % the circles
    N = 10;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rectify images auto on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting image rectification process \n']);tm = clock;
    [rec] = getRectification(imageStack{1},1,0);
    fprintf(['Ending image rectification process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get circles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting circle finding process \n']);tm = clock;
    [circleMatrix MASK tMASK OFFSET] = getCircles(imageStack,rec,N,168);
    fprintf(['Ending circle finding process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map centers to cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting circle mapping process \n']);tm = clock;
    [MAP] = getMap_ver0(tMASK,rec,OFFSET);
    [R] = getBoundingBoxes(MASK);
    for e = 1:numel(R)
        delta = bsxfun(@minus,MAP.centers,R(e).Centroid);
        [J,sidx(e)] = min(sum(delta.*delta,2));
    end
    R(sidx) = R;
    isidx = 1:numel(R);
    isidx = isidx(sidx);
    fprintf(['Ending circle mapping process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    
    para.scales.value = 9;
    para.resize.value = .75;
    BORDER = 100;
    
    clear binarySig rawSignal thresh
    binarySig = zeros(numel(imageStack)-1,numel(R));
    rawSignal = binarySig;
    thresh = zeros(size(numel(R)));
    [pth,nm,ext] = fileparts(imageStack{1});
    fidx = strfind(pth,filesep);
    for e = 1:numel(R)
        tm = clock;
        tmpBB{e} = R(e).BoundingBox;
    end
    [miniStack,miniMask] = diskCrop(imageStack,MASK,tmpBB,BORDER,rec);
    miniStack = single(miniStack);
    miniMask = single(miniMask);
    
    if ~isempty(oPath)
        outFile = [oPath filesep pth((fidx(end)+1):end) '_' num2str(e) '.mat'];
        spoolToDisk_emergence(outFile,miniStack,miniMask);
    end
    fprintf(['Time total estimate:' num2str(numel(R)*etime(clock,tm)/12/60) '\n']);
    close all
end
%{

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) scan for images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataPath = '/iplant/home/jgustin/maizeData/coleoptileEmergence/20170411_Camera4';
    CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
    [o,r] = system(CMD);
    [r] = parseRecords(r);  
    for e = 1:numel(r)
        FileList{e} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        [~,n] = fileparts(FileList{e});
        nm(e) = str2num(n);
    end
    [~,sidx] = sort(nm);
    FileList = FileList(sidx);
    [FileList] = issueBulkTicket(FileList);
    [miniStack,miniMask] = processEmergenceStack_v3(FileList,[]);





    func = cFlow('processEmergenceStack_v3');
    func.setMCRversion('v840');
    func.setMemory(8000);
    [miniStack{1},miniMask{1}] = func(FileList,[]);
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) scan for images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FilePath = '/home/nate/Downloads/Overhead_Compilation/';
    FilePath = '/home/nate/Downloads/Angle_Compilation/';
    FilePath = '/home/nate/Downloads/emergance/';
    FilePath = '/home/nate/Downloads/20151222_Camera1/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/20160106_Camera1/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/Checkerboard/';
    FilePath = '/mnt/scratch1/phytoM/flashProjects/workWithGustin/20160115_Camera1/';
    FilePath = '/mnt/snapper/nate/forEmergance/20160928_Camera3/';
    FilePath = '/mnt/snapper/nate/forEmergance/20161110_Camera3/';
    FilePath = '/mnt/snapper/nate/forEmergance/20161202_Camera3/';
    FilePath = '/home/nate/Downloads/20161219_Camera3/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170131_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170131_Camera4/';

    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera1/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera2/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170220_Camera4/';


    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera1/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera2/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera3/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170317_Camera4/';

    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera1/';
    FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera2/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera3/';
    %FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/20170411_Camera4/';
    FileList = {};
    FileExt = {'tiff','TIF','tif','JPG','jpg'};
    FileExt = {'tiff'};
    verbose = 1;
    FileList = gdig(FilePath,FileList,FileExt,verbose);
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0) viewDrift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = {};
    for e = 1:numel(FileList)
        [p n{e} ex] = fileparts(FileList{e});
        n{e} = str2num(n{e});
    end
    [n sidx] = sort(cell2mat(n));
    FileList = FileList(sidx);

    tmp = imread(FileList{2});
    [J,BOX] = imcrop(tmp);
    ROW = [BOX(2) BOX(2) + BOX(4)];
    COL = [BOX(1) BOX(1) + BOX(3)];

    IJ = [];
    for e = 1:numel(FileList)
        tmp = imread(FileList{e},'PixelRegion',{round(ROW) round(COL)});
        %tmp = imcrop(tmp,BOX);
        IJ(:,:,:,e) =tmp;
        imshow(tmp,[]);
        title(num2str(e));
        drawnow
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) checkerboad check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %checkerBoardChecker(FileList);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) sort the images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    oPath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
    n = {};
    for e = 1:numel(FileList)
        [p n{e} ex] = fileparts(FileList{e});
        n{e} = str2num(n{e});
    end
    [n sidx] = sort(cell2mat(n));
    FileList = FileList(sidx);
    FileList(204:end) = [];
    processEmergenceStack_v2(FileList,oPath);










%}
    %}