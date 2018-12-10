function [tform,resEstimate] = checkerBoardAnalysis(I,box,disp)
    try
        if ischar(I)
            I = double(imread(I))/255;
        end

        
        [tform,resEstimate,imagePoints] = getTform(I);
        

        if disp
            
            imshow(I,[]);
            hold on
            rectangle('Position',box,'EdgeColor','c','LineWidth',3)
            plot(imagePoints(:,1),imagePoints(:,2),'b*')
            plot(imagePoints(:,1),imagePoints(:,2),'ro')
            hold off
        end
    catch ME
        tform = [];
        resEstimate = [];
        ME;
    end

end

function [tform,resEstimate,imagePoints] = getTform(I)
    % hard code checker board size
    cbs = 2.45;
    
    [imagePoints,boardSize] = detectCheckerboardPoints(I);
    
    %{
    imagePoints = cat(3,imagePoints,imagePoints);

    worldPoints = generateCheckerboardPoints(boardSize,cbs);
    imageSize = [size(I,1) size(I,2)];
    params = estimateFisheyeParameters(imagePoints,worldPoints,imageSize,'EstimateAlignment',true);
    I = undistortFisheyeImage(W,params.Intrinsics);

    I = double(I)/255;

    [imagePoints,boardSize] = detectCheckerboardPoints(I);
    %}
    
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