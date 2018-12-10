function [K IN rK] = ker1(IN)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select which D from the cell
    idx = size(IN.D,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % allocate the kernel
    K = zeros(size(IN.D{idx},1),size(IN.D{1},1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calc    
    for i = 1:size(IN.D{idx},1)
        for j = 1:size(IN.D{1},1)
            K(i,j) = mvnpdf(IN.D{1}(j,:),IN.D{idx}(i,:),IN.P);
        end
    end
    if IN.nor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain offsets
        cU = mean(K,2);

        if IN.op == 1
            rU = IN.rU;
            oU = IN.oU;
        else
            rU = mean(K,2);
            oU = mean(K(:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store
        if IN.op == 0
            IN.rU = rU;
            IN.oU = oU;
        end    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calc and normalize
        for i = 1:size(IN.D{idx},1)
            for j = 1:size(IN.D{1},1)
                K(i,j) = K(i,j) - rU(j) - cU(i) + oU;
            end
        end
    end
end

    