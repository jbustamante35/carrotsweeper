function [tf,tform,resEstimate] = isCheckerBoard(oI,tform,resEstimate,GMM,checkerBoardLabelNumber,threshold,disp)
    try
        oI = double(oI);

        sz = size(oI);
        %I = reshape(oI,[prod(sz(1:2)) sz(3)]);
        %I = double(I);
        
        
        label = labelImage(double(oI),GMM);
        %{
        label = GMM.cluster(I);
        
        %}
        label = reshape(label,sz(1:2));
        checkerBoard = label == checkerBoardLabelNumber;



        hcheckBoard = imfill(checkerBoard,'holes');
        ccheckBoard = hcheckBoard == 1  & checkerBoard == 0;
        ccheckBoard = bwlarge(ccheckBoard);
        checkA = sum(ccheckBoard(:));


        tf = checkA > threshold;





        if tf
            [tform,resEstimate,imagePoints] = getTform(oI/255);
        end

        if disp
            out = flattenMaskOverlay(oI/255,checkerBoard,.8,'b');
            out = flattenMaskOverlay(out,ccheckBoard,.2,'b');
            imshow(out,[]);
            if tf
                hold on
                plot(imagePoints(:,1),imagePoints(:,2),'b*')
                plot(imagePoints(:,1),imagePoints(:,2),'ko')
            end
        end
    catch ME
        ME;
    end

end

function [tform,resEstimate,imagePoints] = getTform(I)
    % hard code checker board size
    cbs = 2.45;
    
    [imagePoints,boardSize] = detectCheckerboardPoints(I);

    [imagePoints,resEstimate] = sortPoints(imagePoints,5);

    worldPoints = generateCheckerboardPoints([6 6],cbs);
   
    worldPoints = bsxfun(@minus,worldPoints,mean(worldPoints,1));
    
    [worldPoints,~] = sortPoints(worldPoints,5);

    tform = fitgeotrans(worldPoints,imagePoints,'projective');
            
    resEstimate = resEstimate / cbs;
end

function [pointList,resEstimate] = sortPoints(pointList,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % order the checkerboard points
    kidx = kmeans(pointList(:,2),N);
    L = [];
    for g = 1:N
        L = [L ,kidx==g];
        v(g) = mean(pointList(kidx==g,2));
    end
    [v,sidx] = sort(v);
    L = L(:,sidx);
    resEstimate = mean(diff(v,1,2));
    % my sort
    iS = [];
    for g = 1:N
        fidx = L(:,g)==1;
        d = pointList(fidx,:);
        [~,sidx] = sort(d(:,1));
        d = d(sidx,:);
        iS = [iS;d];
    end
    pointList = iS;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end