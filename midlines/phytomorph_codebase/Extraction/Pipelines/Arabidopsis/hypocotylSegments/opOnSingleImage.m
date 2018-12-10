function [L,WID,dB,P,loc] = opOnSingleImage(imageName,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get mask and skeleton
    I = imread(imageName);
    % read, double, normalize
    tmp = double(I)/255;
    % close the image
    BK = imclose(tmp,strel('disk',61,0));
    % filter and smooth the closed image
    BK = imfilter(BK,fspecial('gaussian',[51 51],9),'replicate');
    % filter the orginal
    tmp = imfilter(tmp,fspecial('gaussian',[51 51],9),'replicate');
    % subtract off the background
    tmp = tmp - BK;
    % normalize
    tmp = bindVec(tmp);
    % threshold value
    level = graythresh(tmp);
    % threshold image
    MASK = double(tmp) < level;
    % remove greater than 1000 objects
    MASK = bwareaopen(MASK,1000);
    % fill holes
    MASK = imfill(MASK,'holes');
    % clear border
    MASK = imclearborder(MASK);
    % get skeleton
    SKEL = bwmorph(MASK,'skel',inf);
    % store skeleton in tmp
    tmp = SKEL;
    % process skeleton 
    for k = 1:20
        ep = bwmorph(tmp,'endpoints');
        tmp(find(ep)) = 0;
    end
    % store skeleton
    SKEL = tmp;
    
    
    
    
    %% get the boundary of the mask
    dB = bwboundaries(MASK);
    for c = 1:numel(dB)
        % get the poly2 mask
        W = poly2mask(dB{c}(:,2), dB{c}(:,1), size(MASK,1),size(MASK,2));
        % select the skeleton based on mask
        W = SKEL.*W;
        % find the skeleton
        [iy ix] = find(W);
        % measure the curvature
        [out] = cwtK_closed_imfilter(dB{c},{[51]});
        % stack the Kurvature
        sig = [out.K;out.K;out.K];
        % dilate the stacked curvature
        peaks = imdilate(-sig,strel('disk',50,0)) == -sig;
        % get the str
        str = size(out.K,1)+1;
        % get the stop
        stp = str + size(out.K) -1;
        % clip out the potential peaks
        peaks = peaks(str:stp);
        % find the peaks
        pidx{c} = find(peaks);
        % get the kurvature values 
        values = out.K(pidx{c});
        % sort the peaks
        [~,sidx] = sort(values);
        pidx{c} = pidx{c}(sidx(1:2));
        loc{c} = dB{c}(pidx{c},:);         
        DP = [fliplr(loc{c})' [ix iy]'];
        T = Radjacency(DP,100);
        [path , pathcost]  = dijkstra(T , 1 , 2);
        P{c} = DP(:,path);
        %fprintf(['Done with analysis: ' num2str(c) ' and image ' num2str(e) '\n']);
    end
    
    
    %% measure
    try
        L = [];
        WID = [];
        tmp = bwdist(~MASK);
        for con = 1:numel(dB)
            IDX = sub2ind(size(tmp),P{con}(2,:),P{con}(1,:));
            values = tmp(IDX);
            WID(con) = mean(values);
            dL = diff(P{con},1,2);
            dL = sum(dL.*dL,1).^.5;
            L(con) = sum(dL);
        end
    catch ME
        ME
    end
    %% display
    if disp
        %close all
        Z = zeros(size(MASK,1),size(MASK,2));
        RGB = cat(3,SKEL,Z,MASK);        
        %imshow(double(RGB),[]);
        image(double(RGB));
        hold on
        for c = 1:numel(dB)
            plot(dB{c}(:,2),dB{c}(:,1),'g')
            plot(loc{c}(:,2),loc{c}(:,1),'r*')
            plot(P{c}(1,:),P{c}(2,:),'y')
        end
        hold off
        drawnow
    end
end