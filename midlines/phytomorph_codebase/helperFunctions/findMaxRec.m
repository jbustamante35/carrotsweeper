function [rec] = findMaxRec(im)
% ----------------------------------------------------------
    close all
    b = imresize(im,[size(im,1) size(im,1)]); %Make the mask square by resizing it by its aspect ratio.
    SC = 1;  % Put 2..4 to scale down the image an speed up the algorithm

    [x1,y1,l1] = findRect(b,SC);          % Lunch the dyn prog algorithm
    [x2,y2,l2] = findRect(rot90(b),SC);   % rotate the image by 90deg and solve
    % Rotate back: x2,y2 according to rot90
    tmp = x2;
    x2 = size(im,1)/SC-y2-l2;
    y2 = tmp;

% Select the best solution of the above (for the original image and for the rotated by 90degrees
        if (l1>=l2)
           corn = sqCorn(x1,y1,l1);
        else
           corn = sqCorn(x2,y2,l2);
        end

        b = imresize(b,1/SC);
%{
    figure;imshow(b>0); hold on;
    plot(corn(1,:),corn(2,:),'O')
%}
    corn = corn*SC;
    corn(1,:) = corn(1,:)*size(im,2)/size(im,1);
%{
    figure;imshow(im); hold on;
    plot(corn(1,:),corn(2,:),'r*')
%}

    rec = [corn(1,1) corn(2,1) corn(1,4)-corn(1,1) corn(2,4)-corn(2,1)];

%rec = 1;
end





function corn = sqCorn(x,y,l)
     corn = [x,y;x,y+l;x+l,y;x+l,y+l]';
end


% ----------------------------------------------------------
function [x,y,l] = findRect(b,SC)
    b = imresize(b,1/SC);
    res = zeros(size(b,1),size(b,2),3);
    % initialize first col
    for i = 1:1:size(b,1)
        if (b(i,1) > 0)
           res(i,1,:) = [i,1,0];
        end
    end
    % initialize first row
    for i = 1:1:size(b,2)
        if (b(1,i) > 0)
           res(1,i,:) = [1,i,0];
        end
    end

    % DynProg
    for i = 2:1:size(b,1)
    for j = 2:1:size(b,2)
        isWhite = b(i,j) > 0;
        if (~isWhite)
           res(i,j,:)=res(i-1,j-1,:); % copy
        else
            if (b(i-1,j-1)>0)  % continuous line
               lineBeg    = [res(i-1,j-1,1),res(i-1,j-1,2)]; 
               lineLenght = res(i-1,j-1,3); 
               if ((b(lineBeg(1),j)>0)&&(b(i,lineBeg(2))>0))  % if second diag is good
                  res(i,j,:) = [lineBeg,lineLenght+1];
               else
                  res(i,j,:)=res(i-1,j-1,:); % copy since line has ended
               end
            else
               res(i,j,:) = [i,j,0];         % Line start
            end
        end

    end
    end

    % check last col
    [maxValCol,WhereCol] = max(res(:,end,3));
    % check last row
    [maxValRow,WhereRow] = max(res(end,:,3));

    % Find max
    x= 0; y = 0; l = 0;
    if (maxValCol>maxValRow)
       y = res(WhereCol,end,1);
       x = res(WhereCol,end,2);
       l = maxValCol;
    else
       y = res(end,WhereRow,1);
       x = res(end,WhereRow,2);
       l = maxValRow;
    end

        corn = [x,y;x,y+l;x+l,y;x+l,y+l]';
    %    figure;imshow(b>0); hold on;
     %   plot(corn(1,:),corn(2,:),'O')
    return;
end



