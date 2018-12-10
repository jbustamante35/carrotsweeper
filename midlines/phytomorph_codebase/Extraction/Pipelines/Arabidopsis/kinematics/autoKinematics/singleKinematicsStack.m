function [F] = singleKinematicsStack(fileList,disp,oPath,DST,DSS)
    fprintf(['****************************************************************************************\n']);
    versionString = ['Starting straight root kinematcs analysis algorithm. \nPublication Version 1.0 - Friday, June 3, 2016. \n'];
    fprintf(versionString);
    fprintf(['****************************************************************************************\n']);
    totalTimeInit = clock;
    initIrods();
    I = imread(fileList{1});
    I = I(:,:,1);
    [midlineM,rootWidth] = getMidline(I);
    [midlineM] = averageMidlines(midlineM);
    [pointList] = trackMidline(midlineM,fileList,rootWidth,DST,DSS,0);
    [domain,vP] = measureVelocityProfile(pointList,1,1,DSS);
    F = cat(3,domain,vP);
    
    if disp
        imshow(I,[]);
        hold on
        plot(path(2,:),path(1,:),'r','LineWidth',3)
    end
    
    if ~isempty(oPath)
        mkdir(oPath);
        [pth,nm,ext] = fileparts(fileList{1});
        nm = fileList{1};
        r = [num2str(round(rand(1)*10^6)) num2str(round(rand(1)*datenum(clock)))];
        save([oPath filesep num2str(r) '_output.mat'],'pointList','F','nm');
    end
    
    close all
    fprintf(['****************************************************************************************\n']);
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
    fprintf(versionString);
    fprintf(['****************************************************************************************\n']);
end

%{

    FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/';
    %FilePath = '/home/nate/Downloads/Wassim''s compilation of old data/';
    %FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/straight/hypoosmotic/';
    FileList = {};
    FileExt = {'tiff','TIF'};
    verbose = 1;
    SET = sdig(FilePath,FileList,FileExt,verbose);
    


    oPath = ['/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/return/straightKinematicsW'];
    singleKinematicsStack(SET{1},0,oPath,1,10);


    

    oPath = ['/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/return/straightKinematics/'];
    for e = 1:numel(SET)
        try
            fprintf(['**********************************************\n']);
            fprintf(['Starting analysis : ' num2str(e) ':' num2str(numel(SET)) '\n']);
            G{e} = singleKinematicsStack(SET{e},0,oPath,1,10);
            fprintf(['Ending analysis : ' num2str(e) ':' num2str(numel(SET)) '\n']);
            fprintf(['**********************************************\n']);
        catch
        end
    end

    for e = 1:numel(G)
        tmp = extractStrain(G{e}(:,:,end),[30 3]);
        mesh(tmp(:,:,2));
        view([0 90]);
        axis([1 size(tmp,2) 1 size(tmp,1)])
        waitforbuttonpress
    end
%}