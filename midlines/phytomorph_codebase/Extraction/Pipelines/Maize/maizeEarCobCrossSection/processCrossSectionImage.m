function [rowN_ENS] = processCrossSectionImage(fileName,oPath,rPath,toSave)

    try
        initIrods();
        mkdir(oPath);
        fileList = {};
        % get file parts
        [pth,nm,ext] = fileparts(fileName);
        % read image
        I = single(imread(fileName))/255;
        G = rgb2gray(I);
        M = G > graythresh(G);
        M = bwareaopen(M,10000);
        M = imclose(M,strel('disk',110));
        R = regionprops(M,'BoundingBox','Centroid','PixelIdxList');
        % for each object in image - process as cob
        for e = 1:numel(R)
            fprintf(['Starting Analysis' num2str(e) ':' num2str(numel(R)) '\n']);
            Mtmp = zeros(size(M));
            Mtmp(R(e).PixelIdxList) = 1;
            tmpB = R(e).BoundingBox - [150 150 -300 -300];
            tmpI = imcrop(I,tmpB);
            [MAJ(e) MIN(e) rowN(e) KERD(e) F(e,:) ROT(:,:,e)] = processCrossSection(double(tmpI),0);
            fprintf(['Ending Analysis' num2str(e) ':' num2str(numel(R)) '\n']);
        end
        %%%%%%%%%
        % average FFT and pick out doubles
        %%%%%%%%%
        uF = mean(F,1);
        
        [sF,sidx] = sort(uF(1:100),'descend');
        sidx = sidx - 1;
        if (sidx(1)/sidx(2) == .5) || (sidx(2)/sidx(1) == 2)
            sidx(1) = sidx(2);
        end
        rowN_ENS = sidx(1);
        rowN_U = mean(rowN);
        
        % 
        h = image(I);
        hold on
        vP = 1.3;
        axis off
        axis equal
        TH = linspace(-pi,pi,2000);
        LineWidth = 1;
        for e = 1:numel(R)
            fprintf(['Starting Display' num2str(e) ':' num2str(numel(R)) '\n']);
            SIG = 50*cos(TH*rowN(e));
            %Ro = .5*(MIN(e)+MAX(e));
            minE(1) = (MAJ(e) - .5*KERD(e));
            minE(2) = (MIN(e) - .5*KERD(e));
            Ro = KERD(e);
            xl = (SIG+minE(1)).*cos(TH);
            yl = (SIG+minE(2)).*sin(TH);
            tmp = (ROT(:,:,e)*[xl(:) yl(:)]')';
            xl = tmp(:,1) + R(e).Centroid(1);
            yl = tmp(:,2) + R(e).Centroid(2);
            plot(xl,yl,'r','LineWidth',LineWidth);

            xl = MAJ(e).*cos(TH);
            yl = MIN(e).*sin(TH);
            tmp = (ROT(:,:,e)*[xl(:) yl(:)]')';
            xl = tmp(:,1) + R(e).Centroid(1);
            yl = tmp(:,2) + R(e).Centroid(2);
            plot(xl,yl,'b','LineWidth',LineWidth);

            xl = (MAJ(e)-KERD(e)).*cos(TH);
            yl = (MIN(e)-KERD(e)).*sin(TH);
            tmp = (ROT(:,:,e)*[xl(:) yl(:)]')';
            xl = tmp(:,1) + R(e).Centroid(1);
            yl = tmp(:,2) + R(e).Centroid(2);
            plot(xl,yl,'b','LineWidth',LineWidth);
            TI = ['Kernel Row Number:' num2str(rowN(e)) '**' 'Kernel Depth:' num2str(KERD(e)) '**' 'Kernel MAJ Width:' num2str(MAJ(e)) '**' 'Kernel MIN Width:' num2str(MIN(e))];
            title(TI);
            %{
            text(R(e).Centroid(1)-100,R(e).Centroid(2),['Kernel Row Number:' num2str(rowN(e))],'Background','w');
            text(R(e).Centroid(1)-100,R(e).Centroid(2)+150,['Kernel Depth:' num2str(MAX(e) - MIN(e))],'Background','w');
            text(R(e).Centroid(1)-100,R(e).Centroid(2)+300,['Kernel Width:' num2str(MAX(e)/rowN(e))],'Background','w');
            %}
            
            
            text(R(e).Centroid(1) - vP*MAJ(e)-100,R(e).Centroid(2),[num2str(e)],'Background','w');
            
            vMIN = R(e).Centroid - vP*MAJ(e);
            vMAX = R(e).Centroid + vP*MAJ(e);
            axis([vMIN(1) vMAX(1) vMIN(2) vMAX(2)]);
            drawnow
            
            % save image file
            fileList{end+1} = [oPath nm '--subNumber--' num2str(e) '.tif'];
            saveas(h,fileList{end});
            fprintf(['Ending Display' num2str(e) ':' num2str(numel(R)) '\n']);
        end
        %close all
        axis([0 size(I,2) -200 size(I,1)]);
        title([]);
        axis equal
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if toSave
            fprintf(['starting save/push phase \n']);
            % save csv data
            fileList{end+1} = [oPath nm '--compile_results.csv'];
            DATA = [MAJ' MIN' KERD' rowN'];
            csvwrite(fileList{end},DATA);
            
            fileList{end+1} = [oPath nm '--row_number_totals.csv'];
            DATA = [rowN_ENS];
            csvwrite(fileList{end},DATA);
            
            
            pushToiRods(rPath,fileList);
            
            
            
            fprintf(['ending save phase \n']);
            
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        close all
    catch ME
        getReport(ME)
        close all
    end
end

%{    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get file list, transform file list for URL, issue tickets for write
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    user = 'garf0012';
    [FileList] = ScanAndIssueNewFilesOniRods(user,'crossSection','maize',{'tif','TIF','tiff','nef'},0);
    numJobs = numel(FileList);    
    remoteOutputLocation = ['/iplant/home/' user '/#plantType#/return/#tissueType#/'];
    remoteOutputLocation = strrep(remoteOutputLocation,'#plantType#','maizeData');
    remoteOutputLocation = strrep(remoteOutputLocation,'#tissueType#','crossSectionData');
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numJobs,'write');


    %{
    fnS = '78010_005'
    fnS = '4N506_005';
    fnS = 'A682_005';
    fnS = '78010_003';
    fnS = 'MN01-150909-cross-0003_mid';
    for e = 1:1000
        [p,nm,ext] = fileparts(FileList{e});
        if strcmp(nm,fnS)
            idxS = e;
        end
    end
    processCrossSectionImage(FileList{idxS},'/mnt/spaldingdata/nate/junkCS/',[],0);
    %}






    func = cFlow('processCrossSectionImage');
    func.setMCRversion('v840');
    func.setMemory(2000);
    for e = 1:numel(FileList)
        fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
        func(FileList{e},'./output/',remoteOutputLocation,1);
        fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    end
    func.submitDag(500,500);

    
    
    processCrossSectionImage(FileList{e},oPath);


    FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/garf0012/crosssectionData/';
    FileList = {};
    FileExt = {'tif'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/garf0012/return/crossSectionData/output/';
    parfor e = 1:numel(FileList)
        processCrossSectionImage(FileList{e},oPath);
    end



    FileList{1} = '/home/nate/Downloads/2016_07_29_sweetness_x_600dpi.png';
    FileList{1} = '/home/nate/Downloads/2016_07_29_trinity_x_600dpi.png';
    oPath = '';
    processCrossSectionImage(FileList{1},oPath);
%}