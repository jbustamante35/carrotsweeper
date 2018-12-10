function [subI] = sampleFAtDomains(I,T,D,SZ,disp,figH)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % I:= image to sample
    % T:= transformation
    % D:= domain - cell array of domains
    % SZ:= domain size for reshape
    %%%%%%%%%%%%%%%
    % subI:= {1} - tranformation for domain
    % subI:= {2} - function sampled at domain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate cell array
    subI = cell(numel(D)+1,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if T is displacement only
    if numel(T) == 2
        dX = T(:)';
        T = eye(3);
        T(1:2,3) = dX;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make D affine if needed
    for d = 1:numel(D)
        if size(D{d},2) ~= size(T,2)
            D{d} = [D{d} ones(size(D{d},1),1)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stack transformation/location for sampling
    subI{1} = T;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each domains
    for d = 1:numel(D)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % transform D to tmpD
        tmpD = (T*D{d}')';
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % clear tmpF
        tmpF = zeros([size(tmpD,1) size(T,3)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % for each channel
        for k = 1:size(I,3)
            % interpolation
            tmpF(:,k) = ba_interp2(I,tmpD(:,2),tmpD(:,1));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % display
        if disp
            figure(figH);
            imshow(I,[]);
            hold on
            plot(tmpD(:,2),tmpD(:,1),'.');
            hold off
            drawnow
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape sample patches
        tmpF = reshape(tmpF,[SZ{d} size(tmpF,2)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % store tensor for freeeze and transport
        subI{1+d} = tmpF;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stack and freeze tensor
    subI = freezeTensor(subI);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end






