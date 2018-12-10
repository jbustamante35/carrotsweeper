function [] = squareStitch(I,L)
    % line - index - [pointIndex,dimIndex,lineIndex,setIndex]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort lines and points
    % first set is horizontal
    % second set is vertical
    % horizontal lines are sorted left to right 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for setIndex = 1:size(L,4)
        
        
        %%%%%%%%%%%%%%%%
        % if horizontal lines then sort points on x
        if setIndex == 1;pointSortDim = 1;else;pointSortDim = 2;end
        % sort the points on a line
        for lineIndex = 1:size(L,3)
            [~,sidx] = sort(L(:,pointSortDim,lineIndex,setIndex),'ascend');
            L(:,:,lineIndex,setIndex) = L(sidx,:,lineIndex,setIndex);
        end
        %%%%%%%%%%%%%%%%
        % if horizontal lines then sort lines on y
        if setIndex == 1;pointSortDim = 2;else;pointSortDim = 1;end
        % sort the lines in a set
        [~,sidx] = sort(L(1,pointSortDim,:,setIndex),'ascend');
        L(:,:,:,setIndex) = L(:,:,sidx,setIndex);
        
        
        
    end
        
    imshow(I,[]);
    hold on
    cnt = 1;
    for point = 1:size(L,1)
        for line = 1:size(L,3)
            for set = 1:size(L,4)
                plot(L(point,1,line,set),L(point,2,line,set),'r.')
                text(L(point,1,line,set),L(point,2,line,set),num2str(cnt))
            end
        end
    end
    
    %lineIntersection(L(1,:,1,
end