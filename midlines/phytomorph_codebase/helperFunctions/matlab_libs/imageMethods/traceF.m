function [G] = traceF(P,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % force matching
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           G  = Point Stack
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           P  = List of coordinates of corners
    %           disp = display (1 = on , 0 = off)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % init the point stack
    G = P{1};
    % for each frame - starting at two
        if disp;
            wbar = waitbar(0,'Matching Points...');
        end
        
    for i = 2:size(P,2)
        % for each point
        for j = 1:size(G,1)
            clear D
            % measure distance to all points in next frame
            for k = 1:size(P{i},1)
                dis = G(j,:,i-1) - P{i}(k,:);
                D(k) = norm(dis);
            end
            [J idx] = min(D);
            G(j,:,i) = P{i}(idx,:);
        end
        fprintf(['Done with matching:' num2str(i) ':' num2str(size(P,2)) '\n'])
            
            if disp;
                wbar = waitbar(i/size(P,2),wbar,'Matching Points...');
            end
    end
    if disp;close(wbar);end
end

