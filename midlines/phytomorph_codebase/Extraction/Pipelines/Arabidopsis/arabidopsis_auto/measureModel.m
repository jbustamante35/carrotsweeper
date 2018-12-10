function [K] = measureModel(U,E,C,sz,frameModel)
    % U :  mean
    % E : frame
    % C : components
    % sz : size of model
    M = PCA_BKPROJ(C,E,U);
    
    X1S = M(:,1:sz(1)*sz(2))';
    X2S = M(:,(sz(1)*sz(2)+1):end)';

    X1S  = reshape(X1S,[sz(1:end-1) 1]);
    X2S  = reshape(X2S,[sz(1:end-1) 1]);
    
    
    compS = cat(1,shiftdim(X1S,-1),shiftdim(X2S,-1));
    
    
    T = diff(compS,1,3);
    lT = sum(T.*T,1).^.5;
    T = bsxfun(@times,T,lT.^-1);
    angle = squeeze(atan2(T(2,:,:,:),T(1,:,:,:)));
    K = diff(angle,1,2);
    
    
    %{

    for t = 1:size(compS,2)
        for e = 1:size(compS,3)
            temp{s}(:,t,e) = frameModel(:,:,t)*[squeeze(compS(:,t,e));1];
        end
    end

    
    
    CL = {'r' 'b' 'g'};
    for t = 1:size(temp{1},2)
        for s = 1:numel(temp)
            midline = squeeze(temp{s}(:,t,:));
            s
            plot(midline(1,:),midline(2,:),CL{s});
            hold all
            quiver(affine{s}(1,3,t),affine{s}(2,3,t),affine{s}(1,1,t),affine{s}(2,1,t),100);
            axis equal
        end
        drawnow
        hold off
    end
    %}
    
end