function [pheno mainPath] = myTraceLow_engine(I,X,toC,toLAT,disp,EXTRA)
    % make the RSA object
    RSA = phytoRSA();
    

    % init the lats to nothing    
    pheno.LAT = [];
    
    % generate box
    BOX = generateCropBox(imfinfo(I),X,EXTRA);
    % read image
    I = imread(I);
    % crop image
    I = imcrop(I,BOX);    
    if size(I,3) == 3;I = rgb2gray(I);end
    if toC;I = imcomplement(I);end
    
    
    % get surface kurvature
    para.scales.value=2;    
    para.resize.value=1;
    K = surKur(I,para);
    
    % threshold the stack on curvature
    tmp = bindVec(K(:,:,1));
    level = graythresh(tmp);
    
    % loop over threshold values if needed
    flag = 1;
    while flag
        % make mask from K
        B = tmp > level;
        % get skeleton
        SKEL = bwmorph(B,'skel',inf);
        % prune
        for k = 1:2
            ep = bwmorph(SKEL,'endpoints');
            SKEL(find(ep)) = 0;
        end
        % get skeleton
        SKEL = bwmorph(SKEL,'skel',inf);
        
        % find skeleton
        [x y] = find(SKEL);
        % construct adjacency
        DP = [x y]';
        T = Radjacency(DP,3);
        
        % recenter the clicks
        x1 = X(:,1) - BOX(1);
        x2 = X(:,2) - BOX(2);

        % single trace    
        [idx(1)] = snapTo(DP',[x2(1) x1(1)]);
        [idx(2)] = snapTo(DP',[x2(2) x1(2)]);
        [path , pathcost]  = dijkstra(T , idx(1) , idx(2));
        pheno.gamma = DP(:,path);
        pheno.BOX = BOX;
        
        
        mainPath = phytoPath(RSA,DP(:,path));
        RSA.addBranch(mainPath);
        
        flag = 0;
        if isempty(path)
            level = level - .05*level;
            flag = 1;
        end
        
        
        
        
        if toLAT
            try
                % make root map
                midx = sub2ind(size(B),round(pheno.gamma(1,:)),round(pheno.gamma(2,:)));
                rM = zeros(size(B));
                rM(midx) = 1;
                rM = imdilate(rM,strel('disk',2,0));

                % construct branch map
                BM = imfilter(double(SKEL),ones(3),'replicate');
                BM = BM.*SKEL.*rM;
                BM = BM>=4;
                R = regionprops(BM,'Centroid');
                BMi = [];
                for e = 1:numel(R)
                    BMi(e) = snapTo(pheno.gamma',fliplr(R(e).Centroid));
                end
                BMi = path(BMi)';
                BMi = [BMi idx];

                % construct tip map
                TM = bwmorph(SKEL,'endpoints');
                [tipPoints1 tipPoints2] = find(TM);
                for e = 1:numel(tipPoints1)
                    TMi(e) = snapTo(DP',[tipPoints1(e) tipPoints2(e)]);
                end

                % trace lats
                LAT = traceLats(DP,T,BMi,TMi,inf);
                
                %  trim lats - if any lat roots snap to the ends of the
                %  main root - then remove
                LAT = trimLats(LAT,idx);
                
                
                for e = 1:numel(LAT)
                    child = phytoPath(RSA,LAT(e).gamma);
                    mainPath.attachChild(child,LAT(e).bpIndex);
                end                
                % give the lats no lats
                for e = 1:numel(LAT)
                   LAT(e).LAT = [];
                end
                pheno.LAT = LAT;
            
            catch ME
                   ME 
            end
            
            end
                
        
    end
    
    
    % quick display
    if disp
        close all
        imshow(I,[])
        imshow(cat(3,SKEL,B,tmp),[]);
        hold on
        for e = 1:numel(pheno)
            plot(pheno(e).gamma(2,:),pheno(e).gamma(1,:),'r')
        end
        if toLAT
            plot(DP(2,BMi),DP(1,BMi),'r*');
            plot(tipPoints2,tipPoints1,'b*');
            for e = 1:numel(LAT)
                plot(LAT(e).gamma(2,:),LAT(e).gamma(1,:),'k')
            end
        end
        drawnow
    end
    
    T = [eye(2) flipud(BOX(1:2)')];
    RSA.moveRSA(T);
    
    % recenter the root to global frame
    pheno.gamma = bsxfun(@plus,pheno.gamma,flipud(pheno.BOX(1:2)'));
    % recenter lats too
    if toLAT
        for l = 1:numel(pheno.LAT)   
            pheno.LAT(l).gamma = bsxfun(@plus,pheno.LAT(l).gamma,flipud(pheno.BOX(1:2)'));
        end
    end


end

function [BOX] = generateCropBox(info,X,EXTRA)
    UL = X(2,:) - EXTRA;
    if UL(1) <= 0
        UL(1) = 1;
    end
    if UL(2) <= 0
        UL(2) = 1;
    end    
    HW = X(1,:) - X(2,:);
    HW = HW + 2*EXTRA;
    BOX = [UL HW];
    if BOX(3) > info.Width
        BOX(3) = info.Width
    end
    if BOX(4) > info.Height
        BOX(4) = info.Height
    end
    
end

function [idx] = snapTo(X,P)
    idx = [];
    delta = bsxfun(@minus,X,P);
    delta = sum(delta.*delta,2);
    [~,idx] = min(delta);
end

function [LT] = traceLats(DP,T,BI,TI,minP)
    cnt = 1;    
    % for each end point - try to connect it to each branch point
    for t = 1:numel(TI)
        path = [];
        pathcost = [];
        % for branch point
        for b = 1:numel(BI)
            [path{b} , pathcost(b)]  = dijkstra(T , TI(t) , BI(b));
        end
        [J,sidx] = min(pathcost);        
        % if the cost is less than a threshold and greater than 0
        if J < minP & J > 0
            LT(cnt).gamma = DP(:,path{sidx});
            LT(cnt).idx = path{sidx};
            LT(cnt).bpIndex = BI(sidx);
            cnt = cnt + 1;
        end
    end
end

function [LT] = trimLats(LT,capsIDX)
    % create default not trim vector
    rm = zeros(numel(LT,1));
    % for each lat root
    for e = 1:numel(LT)
        if any(LT(e).idx == capsIDX(1)) || any(LT(e).idx == capsIDX(2))
            rm(e) = 1;
        end
    end
    LT(logical(rm)) = [];
end
