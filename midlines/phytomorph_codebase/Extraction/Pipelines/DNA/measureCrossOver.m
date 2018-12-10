function [] = measureCrossOver(fileName,oPath,rPath)
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get file name information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract needed parts from file name
        [pth,tmpNM,ext] = fileparts(fileName);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get crop boxes and and histograms
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get chromosome boxes
        [I,BOXES,rM,R] = getChromosomeCropBoxes(fileName);
        % get color histograms
        for k = 1:size(I,3)
            HISTO(:,k) = imhist(I(:,:,k));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display the crop boxes in white
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mainH = figure;
        mainA = axes(mainH);
        image(mainA,I);
        hold on
        % render crop boxes in white on the large red objects
        renderCropBoxOnImage(mainH,BOXES,'w');
    
        
        
        
        h2 = figure;
        a1 = axis;
        mkdir(oPath);

        fileList = {};
        
       

       
        axis off;
        hold on
        drawnow
       
        
        
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over each crop box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:numel(BOXES)
            try
                
                skeletonPoints = [];
                vec = [];

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % crop out the color image and get the first pass at red mask
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % crop out the color image
                tmpI = imcrop(I,BOXES{e});
                % get the mask of the object under inspection
                tmpREDMASK = imcrop(rM,BOXES{e});
                tmpREDMASK = bwlarge(tmpREDMASK);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % crop out the color image and get the first pass at red mask
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                avergeRED = tmpI(:,:,1).*tmpREDMASK;
                avergeRED = sum(avergeRED(:))/sum(tmpREDMASK(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the contour and the skeleton for the largest red object
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                scaleFactor = 4;
                padValue = [10 10];
                meltIter = 50;
                alpha = .75;
                [contourPoints,skeletonPoints] = getContourSkeleton(tmpREDMASK,padValue,scaleFactor,meltIter,alpha);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the end points of the skelton
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [skeletonIDX,doubleMeasure(e,1)] = findChromosomeEndPoints(contourPoints,skeletonPoints);

                doubleMeasure(e,2) = sum(tmpREDMASK(:));


                
                
                
                if ~(doubleMeasure(e,1) > 2*10^5 && doubleMeasure(e,2) > 1200)



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % dijkstra trace the midline
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    T = Radjacency(skeletonPoints',.5);
                    [path , pathcost]  = dijkstra(T , skeletonIDX(1) , skeletonIDX(2));
                    path = skeletonPoints(path,:);
                    %{
                    % display after contour and midline and dijkstra
                    figure;
                    imshow(tmpI,[]);
                    hold on
                    plot(dC(:,1),dC(:,2),'y')
                    plot(skeletonPoints(:,1),skeletonPoints(:,2),'y.')
                    plot(path(:,1),path(:,2),'r')
                    hold off
                    %waitforbuttonpress
                    %}
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % dijkstra trace the midline
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % threshold the blue
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % threshold the blue channel
                    bM = tmpB > graythresh(tmpB);
                    out = flattenMaskOverlay(tmpI,logical(bM),1,'b');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % threshold the blue
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %}




                    [DomainS,DomainG] = extendChromosomeMidline(path,fliplr(BOXES{e}(1:2)),tmpI);
                    dsz = size(DomainG);
                    %{
                    % check domainG
                    imshow(I,[]);
                    hold on
                    for e = 1:size(DomainG,2)
                        plot(DomainG(:,e,2),DomainG(:,e,1),'r')
                    end
                    for e = 1:size(DomainG,1)
                        plot(squeeze(DomainG(e,:,2)),squeeze(DomainG(e,:,1)),'r')
                    end
                    %}



                    vec = [];
                    for k = 1:size(I,3)
                        vec(:,k) = ba_interp2(I(:,:,k),DomainS(:,2),DomainS(:,1));
                    end
                    vec = reshape(vec,[dsz(1) dsz(2) 3]);
                    %waitforbuttonpress
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % sub-sample
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % remake curve co-ordinates based on red mask
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    vecR = vec(:,:,1);
                    vecRM = vecR > graythresh(vecR);
                    vecRM = bwlarge(vecRM);
                    fidx = find(any(vecRM,2));
                    oldN = size(path,1);
                    path = squeeze(DomainG(fidx(1):fidx(end),(end-1)/2,:));
                    addedP = round((size(DomainG,1) - size(path,1))/2);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % remake curve co-ordinates based on red mask
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % flip with blue on top
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % sum blue hori
                    sigB = sum(vec(:,:,3),2);
                    sigBM = zeros(size(vec,1),size(vec,2));

                    % find blue peak
                    [~,idxB] = max(sigB);
                    sigBM((idxB-8):(idxB+8),:) = 1;
                    %hold on
                    %plot(sigB,1:size(vec,1),'b');
                    % sum red vert
                    sigR = sum(vec(:,:,1),1);
                    % make red mask
                    maskStrip = vec(:,:,1) > graythresh(vec(:,:,1));
                    maskStrip = bwlarge(maskStrip);
                    %
                    sigM = sum(maskStrip,2);
                    % make green mask
                    ridx = find(maskStrip);
                    greenChannel = vec(:,:,2);
                    dotMask = greenChannel > graythresh(greenChannel(ridx));
                    greenChannelFilter = maskStrip.*greenChannel.*dotMask;
                    igreenChannelFilter = sum(greenChannelFilter,2);
                    dotMaskPeaks = (imdilate(igreenChannelFilter,strel('disk',11)) == igreenChannelFilter).*igreenChannelFilter;
                    dotMaskPeaks = dotMaskPeaks.*any(maskStrip,2);
                    %dotMask = imdilate(dotMask,strel('disk',5,0));
                    %dotMask = bwareaopen(dotMask,10);
                    gR = regionprops(logical(dotMaskPeaks),'Centroid');
                    
                    % save for learning
                    fileList{end+1} = [oPath tmpNM '_'  num2str(e) '_modelData.mat'];
                    save(fileList{end},'vec','maskStrip');
                    %{
                    rmIDX = [];
                    for g = 1:numel(gR)
                        if gR(g).Centroid(2) > (addedP) & gR(g).Centroid(2) < (addedP + size(path,1))
                            rmIDX(g) = false;
                        else
                            rmIDX(g) = true;
                        end
                    end
                    %}
                    %gR(logical(rmIDX)) = [];
                    %plot(1:size(vec,2),sigR,'r');
                    %hold off
                    %waitforbuttonpress
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % flip with blue on top
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






                    % measure length
                    dP = diff(path,1,1);
                    dL = sum(sum(dP.*dP,2).^.5,2);
                    dL_c = cumsum(sum(dP.*dP,2).^.5,1);


                    if idxB > size(vec,1)/2
                        vec = flipdim(vec,1);
                        idxB = -(idxB - size(vec,1)/2) + size(vec,1)/2;
                        for g = 1:numel(gR)
                            gR(g).Centroid(2) = -(gR(g).Centroid(2)-size(vec,1)/2) + size(vec,1)/2;
                        end
                        gR = flipdim(gR,1);
                    end



                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % make model
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    midx = find(sigM);
                    WID = 50;
                    Model = zeros(size(vec)+[0 100 0]);
                    Model(find(sigM),(end/2-5):(end/2+5),1) = 1;
                    Model((idxB-3):(idxB+3),(end/2-8):(end/2+8),3) = 1;
                    for g = 1:numel(gR)
                        Model((round(gR(g).Centroid(2))-1):(round(gR(g).Centroid(2))+1),(end/2-5):(end/2+5),2) = 1;
                        Model((round(gR(g).Centroid(2))-1):(round(gR(g).Centroid(2))+1),(end/2-5):(end/2+5),1) = 0;
                    end


                    Model((midx(1)-3):(midx(1)-1),(end/2-15):(end/2+15),1) = 1;
                    Model((midx(1)-3):(midx(1)-1),(end/2-15):(end/2+15),2) = 1;
                    Model((midx(end)+1):(midx(end)+3),(end/2-15):(end/2+15),1) = 1;
                    Model((midx(end)+1):(midx(end)+3),(end/2-15):(end/2+15),2) = 1;


                    %Model = padarray(Model,[0 50],0,'both');
                    figure(h2);
                    %imshow([vec Model],[]);
                    image([vec Model]);



                    text(size(vec,2)+size(Model,2)/2-4-30,midx(1)-2,'0','Color','y');
                    text(size(vec,2)+size(Model,2)/2-4-30,midx(end)+2,num2str(round(sum(dL))),'Color','y');


                    text(size(vec,2)+size(Model,2)/2-4-60,idxB,[num2str(round(dL_c(idxB-addedP+1))) '--------------' num2str(dL_c(idxB-addedP+1)/sum(dL),2)],'Color','b');
                    for g = 1:numel(gR)
                        text(size(vec,2)+size(vec,2)-40-WID/2+50,round(gR(g).Centroid(2)),[num2str(round(dL_c(round(gR(g).Centroid(2)-addedP+1)))) '-----' num2str((dL_c((gR(g).Centroid(2)-addedP+1)))/sum(dL),2)],'Color','g');
                    end


                    outputData = {};
                    outputData{1,1} = 'fileName';
                    outputData{2,1} = tmpNM;
                    outputData{1,2} = 'chromosomeLength';
                    outputData{2,2} = round(sum(dL));
                    outputData{1,3} = 'centromere_ABS_Position';
                    outputData{2,3} = (dL_c(idxB-addedP+1));
                    outputData{1,4} = 'centromere_PER_Position';
                    outputData{2,4} = dL_c(idxB-addedP+1)/sum(dL);
                    outputData{1,5} = 'avergeRED';
                    outputData{2,5} = avergeRED;
                    outputData{1,6} = 'numberCrossOvers';
                    outputData{2,6} = numel(gR);
                    outputData{1,7} = 'boxNumber';
                    outputData{2,7} = e;
                    for g = 1:numel(gR)
                        outputData{1,end+1} = ['foci_ABS_Position_' num2str(g)];
                        outputData{2,end} = round(dL_c(round(gR(g).Centroid(2)-addedP+1)));
                        outputData{1,end+1} = ['foci_PER_Position_' num2str(g)];
                        outputData{2,end} = (dL_c((gR(g).Centroid(2)-addedP+1)))/sum(dL);
                        outputData{1,end+1} = ['foci_intensity_' num2str(g)];
                        outputData{2,end} = igreenChannelFilter(gR(g).Centroid(2));
                    end
                    fileList{end+1} = [oPath tmpNM '_'  num2str(e) '_data.csv'];
                    cell2csv(fileList{end},outputData);





                    drawnow
                    fileList{end+1} = [oPath tmpNM '_'  num2str(e) '_straight.tif'];
                    imwrite([vec Model],fileList{end});
                    fileList{end+1} = [oPath tmpNM '_'  num2str(e) '_straight2.tif'];
                    saveas(h2,fileList{end});
                    pause(.001)
                    hold off
                    %waitforbuttonpress



                    figure(mainH)
                    hold on

                    plot(R(e).Centroid(1),R(e).Centroid(2),'y*')
                    rectangle('Position',R(e).BoundingBox,'EdgeColor','r')
                    text(R(e).Centroid(1)+10,R(e).Centroid(2)+10,num2str(e),'Background','w')
                    drawnow
                    % save local image for push to remote
                    fileList{end+1} = [oPath tmpNM  '_whole.tif'];
                    saveas(mainH,fileList{end});
                    % save local mat file - include HISTO
                    fileList{end+1} = [oPath tmpNM  '_meta.mat'];
                    save(fileList{end},'HISTO');
                    %waitforbuttonpress
                    
                    
                    
                end
            catch ME
                getReport(ME)

                %break
            end
                
                
                
          

        end
        figure(mainH);
        hold off
        figure(h2);
        hold off

        %%%%%%%%%%%%%%%%%%%%%%%%%
        % push to iRODS
        %%%%%%%%%%%%%%%%%%%%%%%%%
        pushToiRods(rPath,fileList);
        fprintf(['ending save phase \n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%
         close all 
    catch ME
        getReport(ME)
    end
    
end






%{



    F = readtext('/home/nate/Downloads/5BackwardsImages.csv');
    mkdir('./FLIPoutput/');
    for f = 2:size(F,1)
        fileName = [F{f,3} filesep F{f,1} '.tif'];
        measureCrossOver(fileName,'./FLIPoutput/',[]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run on condor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%
    % dig for data
    %%%%%%%%%%%%%%%%%%%%%%%%
    dataPath = ['/iplant/home/petersoapes/for tests%'];
    dataPath = ['/iplant/home/petersoapes/Peromyscus%'];
    %dataPath = ['/iplant/home/petersoapes/NewMusData%'];
    %dataPath = ['/iplant/home/petersoapes/OLDImageSet%'];
    
    CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
    [o,r] = system(CMD);
    [r] = parseRecords(r);
    FileExt = 'rev.tif';
    FileList = {};
    for e = 1:numel(r)
        [p,nm,ext] = fileparts(r(e).DATA_NAME);
        %if any(strcmp(ext(2:end),FileExt))
        if r(e).DATA_NAME((end-numel(FileExt)):end)
            FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%
    % issue write over return folder
    %%%%%%%%%%%%%%%%%%%%%%%%
    numJobs = numel(FileList);
    remoteOutputLocation = '/iplant/home/petersoapes/returns/returnPaper_Peromyscus/';
    %remoteOutputLocation = '/iplant/home/petersoapes/returns/returnPaper_NEW/';
    %remoteOutputLocation = '/iplant/home/petersoapes/returns/returnPaper_OLD/';
    
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numJobs,'write');
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % issue read tickets over images
    %%%%%%%%%%%%%%%%%%%%%%%%
    [FileList] = issueBulkTicket(FileList);




    
    %% run local
    %measureCrossOver(FileList{e},'./output/',[]);



    %%%%%%%%%%%%%%%%%%%%%%%%
    % launch
    %%%%%%%%%%%%%%%%%%%%%%%%
    func = cFlow('measureCrossOver');
    func.setMCRversion('v930');
    numJobs = 10;
    numJobs = numel(FileList);
    for e = 1:numJobs
        fprintf(['start generating job:' num2str(e) ':' num2str(numJobs) '\n']);
        %mecka(algorithmFlag,FileList{e},numberOfObjects,'./output/',remoteOutputLocation,1,1,imageRES,1);
        func(FileList{e},'./output/',remoteOutputLocation);
        fprintf(['end generating job:' num2str(e) ':' num2str(numJobs) '\n']);
    end


    %{
    if numJobs ~= 0
        [CMD] = uncLog('ph:l:',rawFileList(1:numJobs),'add',['maize' '-' analysisType],'1.0',{'0'},{'1'},1);
    end
    %}
    
    
    
    
    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    auth = strrep(auth,'data.iplantc.org','davos.cyverse.org');


    func.submitDag(auth,150,150);
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run on condor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%}









