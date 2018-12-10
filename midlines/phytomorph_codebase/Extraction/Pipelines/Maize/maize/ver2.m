modelGEN = 1;
for nTAR = 1:27
    %% find files and generate rand number
    pth = '/mnt/scratch5/takeshi_dev_run2/mat/'
    pth = '/mnt/scratch5/takeshi_dev_run_2013.06.15/';
    pth = '/mnt/scratch5/takeshi_dev_run_2013.06.21/';
    pth = '/mnt/spaldingdata/nate/devRuns/takeshi_dev_run_2013.07.03/';
    pth = '/mnt/spaldingdata/nate/devRuns/takeshi_dev_run_2013.08.11/';
    pth = '/mnt/spaldingdata/nate/devRuns/takeshi_dev_run_2013.13.11/';
    pth = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/allExtraction/';
    dataFiles = gdig(pth,{},{'mat'},1);
    numel(dataFiles)
%%  
h = figure;
[D] = modelGeneration(dataFiles,h);
    [target] = expectedKernelNumber(dataFiles{1});
    
    %% gather expected number of kernels from the file name
    for i = 1:numel(dataFiles)
        % find parts of file(s)
        [p n e] = fileparts(dataFiles{i});
        % remove the name
        n(end-4:end) = [];
        % find the key char
        fidx = strfind(n,'_');
        %
        if ~isempty(fidx) & numel(fidx) >= 2
            SNIP = n(1:fidx(2));
            n(1:fidx(2)) = [];
            n = ['_' n '_'];
            fidx = strfind(n,'_');
            target(i) = numel(fidx) - 1;
        else
            target(i) = 0;
        end
    end
    %% loop load
    clear dataS
    % if nTAR == 0 then build model
    if nTAR == 0
        nTAR = 6;
    end
    % find the ones with 6 as the target
    fidx = find(target==nTAR);
    totL = numel(fidx);
    % load the data
    parfor i = 1:totL
        try
            % get number to load
            cnt = fidx(i);
            % start clock
            tic
            % load
            dataS{i} = genLoader(dataFiles{cnt},cnt);        
            % report
            fprintf(['Done loading' num2str(i) ':' num2str(totL) ':' num2str(toc) '\n']);
        catch ME

        end
    end
    % render D,A,L,
    D = [];
    A = [];
    L = [];
    groups = [];
    target = [];
    for i = 1:numel(dataS)
        D = [D;dataS{i}.V];
        A = [A;dataS{i}.A];
        L = [L;dataS{i}.L];
        groups = [groups;dataS{i}.idx];
        target(i) = nTAR;
        imageStack{i} = dataS{i}.imageStack;
    end
    %% normalize data
    for e = 1:size(D,1)
        D(e,:) = bindVec(D(e,:));
    end
    %% model generation
    if modelGEN
        h = figure;
        modelGEN = 0;
        UQg = unique(groups);
        selV = [];
        for i = 1:numel(dataS)

            fidx = find(groups == UQg(i));

            S = dataS{i}.imageStack;
            I = imread(S{1}.fileName);
            SNIP = 400;
            constM = zeros(size(I));

            sz = size(I);
            clear V V1 V2 I J G;
            parfor j = 1:61
                I = imread(S{j}.fileName);        
                V1{j} = mean(I,1);
                V2{j} = mean(I,2);
                J{j} = I;
                [g1 g2] = gradient(double(I));
                g = (g1.^2 + g2.^2).^.5;
                G{j} = g;
            end

            J = cell2mat(J);
            J = reshape(J,[sz 61]);

            G = cell2mat(G);
            G = reshape(G,[sz 61]);

            V1 = cell2mat(V1');
            V2 = cell2mat(V2);

            sV2 = std(V2,1,2);
            sV1 = std(V1,1,1);

            uG = mean(G,3);
            sG = std(G,1,3);

            sJ = std(double(J),1,3);

            sJ = sJ.*sG;

            sJ = normalize(sJ);
            thresh = graythresh(sJ);
            MASK = sJ > thresh;
            MASK = imclose(MASK,strel('disk',21));
            MASK = imerode(MASK,ones(5));

            MASK = bwareaopen(MASK,100);
            CC = bwconncomp(MASK,8);
            R = regionprops(MASK,'Centroid');


            clickPos = [R.Centroid];
            clickPos = reshape(clickPos,[2 numel(R)])';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % snap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sidx = [];
            for k = 1:size(clickPos,1)
                dist = [];
                for l = 1:size(dataS{i}.gamma,1)
                    delta = clickPos(k,:) - dataS{i}.gamma(l,:);
                    dist(l) = norm(delta);
                end
                [J sidx(k)] = min(dist);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if numel(sidx) == 6
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for view points
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % draw data for user
                close all
                I = myReader(dataS{i}.imageStack{1}.fileName,'toGray',1);                
                imshow(I);
                hold on;
                nProps.Color = 'g';
                dataS{i}.gD.pointList{1}.view(h,nProps);   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % view graphs
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for u = 1:numel(sidx)
                    vProps.Color = 'r';
                    dataS{i}.G.sG{sidx(u)}.N{1}{1}.view(h,vProps);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get patches from seletions
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                selV = [selV;D(fidx(sidx),:)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                title(num2str(i));
                drawnow
            end
        end
        
        % after first pass create model
        uM = mean(selV,1);
        uM = bindVec(uM);
        delta = bsxfun(@minus,selV,uM);
        dist = sum(delta.*delta,2).^.5;
        
        % find model index and generate model
        ridx = dist > 6;
        selV(ridx,:) = [];
        
        % after first pass create model
        uM = mean(selV,1);
        %uM = bindVec(uM);
        delta = bsxfun(@minus,selV,uM);
        dist = sum(delta.*delta,2).^.5;
        fidx = dist < 5.5;
        
        [S C U BV L ERR LAM] = PCA_FIT_FULL(selV(fidx,:),50);
        LAM = diag(LAM).^-1;
        LAM = diag(LAM);
    end
    %% mass load switch view
    UQg = unique(groups);
    selV = [];
    VIEW = 1;
    SPOOL = 1;
    LEVEL = 2;
    %outPort.csvPath = '/mnt/scratch5/takeshi_dev_run6/csvSPOOL/';
    outPort.csvPath = '/mnt/spaldingdata/nate/devRuns/takeshi_dev_run_2013.08.11_reduction/csvSPOOL/';
    outPort.csvPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/allExtraction/csvSPOOL/';
    %outPort.imagePath = '/mnt/scratch5/takeshi_dev_run6/imageSPOOL/';
    outPort.imagePath = '/mnt/spaldingdata/nate/devRuns/takeshi_dev_run_2013.08.11_reduction/imageSPOOL/';
    outPort.imagePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/allExtraction/imageSPOOL/';
    mkdir(outPort.csvPath);
    mkdir(outPort.imagePath);
    for i = 1:numel(dataS)
        % get the ith image stack
        S = dataS{i}.imageStack;
        % find the values fro the ith group
        fidx = find(groups == UQg(i));
        % read the first image
        I = imread(S{1}.fileName);  
        %%%%%%%%%%%%%%%%%%%%%%%%
        % get the data from the indices
        subD = D(fidx,:);
        %%%%%%%%%%%%%%%%%%%%%%%%
        % get the data from the indices
        subD1 = reshape(D(fidx,:),[numel(fidx) 20 200]);
        % circshift
        subD1 = circshift(subD1,[0 0 100]);
        subD1 = reshape(subD1,size(subD));
        %%%%%%%%%%%%%%%%%%%%%%%%
        % project the data into the comp space
        C = PCA_REPROJ(subD,BV,U);
        C1 = PCA_REPROJ(subD1,BV,U);
        % sim the data
        SIM = PCA_BKPROJ(C,BV,U);
        SIM1 = PCA_BKPROJ(C1,BV,U);
        % measure the error
        err = SIM-subD;
        err1 = SIM1-subD1;
        % norm of error
        edist = sum(err.*err,2).^.5;
        edist1 = sum(err1.*err1,2).^.5;
        % distance in model space
        dist = sum(C.*C,2);
        dist1 = sum(C1.*C1,2);
        % total distance in normal and tangent space
        tdist = dist + edist;
        tdist1 = dist1 + edist1;
        % cat data
        TOT = [tdist tdist1];
        tdist = min(TOT,[],2);
        % sort
        [sD sidx] = sort(tdist);
        sidx = sidx(1:nTAR);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if view
        if VIEW
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for view points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % draw data for user
            close all
            I = myReader(dataS{i}.imageStack{1}.fileName,'toGray',1);
            h = imshow(I);
            hold on;
            nProps.Color = 'g';
            dataS{i}.gD.pointList{1}.view(h,nProps);   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % view graphs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for u = 1:numel(sidx)
                vProps.Color = 'r';
                dataS{i}.G.sG{sidx(u)}.N{1}{1}.view(h,vProps);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            title(num2str(i));
            drawnow
            pause(.5);
        end


        if SPOOL
            % spool angles
            angle = -A(fidx,:);
            [p n ext] = fileparts(S{1}.fileName);
            pidx = strfind(p,'/');
            p = p(pidx(end-LEVEL)+1:end);
            fn = [outPort.csvPath strrep(p,'/','--') '.csv'];
            [JUNK oidx] = sort(dataS{i}.gamma(sidx,2));
            fn
            csvwrite(fn,angle(sidx(oidx),:));

            % spool images
            fn = [outPort.imagePath strrep(p,'/','--') '.tif'];
            fn
            saveas(h,fn);
        end
        i
    end
end