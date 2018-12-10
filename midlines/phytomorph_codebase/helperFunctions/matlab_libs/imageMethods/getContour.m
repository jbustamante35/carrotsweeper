function [contourFamily] = getContour(BS) 
    try
        cnt = 1;
        % for each binary image find the contour
        for i = 1:size(BS,3)
            % get the binary image and trace the contour
            B = BS(:,:,i);
            x1 = find(B(:,1)); 
            contour = bwtraceboundary(B,[x1(1) 1],'NE',8);
            %fidx = find(contour(:,2) == 1);
            %contour(fidx,:) = [];
            %contour = [[contour(1,1) 1];contour;[contour(end,1) 1]];
            % filter contour on lenght
            if size(contour,1) > 200
                % record the contour for the ith image
                contourFamily(cnt).contour = contour;
                cnt = cnt + 1;
            end
        end
    catch ME
         stop = 1;
    end
       
end