function [sweepT,integrationSweep,petNUM,petDIA,offset] = topSweep(top,saveDisplay,oName,offset)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize sweeper
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % locations for petiole number spots to query for data
    pcLOC = [60 80 100 120 140];
    % height of the grid
    MAX_GRID_LENGTH = 1800;
    % location to start the grid
    START_GRID = 100;
    % number of marks on the display grid along the radius
    DISPLAY_GRID_LENGTH_NUM = 10;
    % number of marks on the display grid along angle
    DISPLAY_GRID_ANGLE_NUM = 10;
    % strech ratio of grid
    GRID_RATIO = 3150/1800;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize sweeper
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make grid for top sweep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [d1,d2] = ndgrid(linspace(0,-pi,1000),linspace(START_GRID,MAX_GRID_LENGTH,MAX_GRID_LENGTH-START_GRID));
    X1 = d2.*cos(d1);
    X2 = GRID_RATIO*d2.*sin(d1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make grid for top sweep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make display grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [di1,di2] = ndgrid(linspace(0,-pi,DISPLAY_GRID_ANGLE_NUM),linspace(START_GRID,MAX_GRID_LENGTH,DISPLAY_GRID_LENGTH_NUM));
    Xi1 = di2.*cos(di1);
    Xi2 = GRID_RATIO*di2.*sin(di1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make display grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find arc length along grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dL = (diff(X1,1,1).^2 + diff(X2,1,1).^2).^.5;
    nL = mean(dL,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if nargin <= 3
        offset = [size(top,2)/2 size(top,1)];
        sig = mean(top(end-20:end,:),1);
        fidx = find(sig > .2);
        offset(1) = mean(fidx);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the sweep - offset is to push the sweep grid to the proper location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sweepT = [];
    for i = 1:size(X1,2)
        sweepT(i,:) = ba_interp2(top,X1(:,i)+offset(1),X2(:,i)+offset(2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % look at select paths along the sweep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the area measure - surrogate for petiole diameter at a sweep location
    sweepTPD = [];
    % number of objects found along the sweep path
    sweepTPN = [];
    % for each defined sweep path
    for i = 1:numel(pcLOC)
        % get the selected sweep path
        sig = sweepT(pcLOC(i),:);      % sample the sweep at petiole locations
        % remove noirs
        sig = bwareaopen(sig,5);
        % threshold and get area
        R = regionprops(sig > .2,'Area');
        % number of "petioles" found at this sweep location
        sweepTPN(i) = numel(R);
        % get the area of each "petiole"
        sweepTPD(i) = mean([R.Area]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make measurements along the grid
    % find the max number of petioles to make the call on total number
    [sweepTPN midx] = max(sweepTPN);
    % another count is the mean of the count - FIX THIS
    petNUM = mean(sweepTPN);
    % measure the diameter as the averge at the selected ring
    petDIA = sweepTPD(midx);
    % integrate along the sweep greid
    integrationSweep = sum(sweepT,2)';
    % flip
    integrationSweep = fliplr(integrationSweep);
    % normalize by the acrlength
    integrationSweep = integrationSweep.*fliplr(nL);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE DISPLAY IMAGES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if saveDisplay
        % make the top image
        image(255*top);
        % set the colormap
        colormap(gray);
        hold on
        % render the display grid
        for i = 1:size(Xi1,2)
            plot(Xi1(:,i)+offset(1),Xi2(:,i)+offset(2),'Color',[1 .5 0],'LineWidth',3);
        end
        % plot select rings
        for i = 1:numel(pcLOC)
            plot(X1(:,pcLOC(i))+offset(1),X2(:,pcLOC(i))+offset(2),'c','LineWidth',2);
        end
        % render the display grid in the ortho-normal direction
        for i = 1:size(Xi1,2)
            plot(Xi1(i,:)+offset(1),Xi2(i,:)+offset(2),'Color',[1 .5 0],'LineWidth',3);
        end
        % title the image
        title(num2str(petNUM(end)));
        plot(d2(1,:)+offset(1),size(top,1)*ones(size(integrationSweep)) - fliplr(integrationSweep),'c','LineWidth',2);
        % format the image
        axis off
        axis equal
        % save the image
        saveas(gca,oName);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE DISPLAY IMAGES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    close all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end