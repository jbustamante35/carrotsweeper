function [IN] = MCA(IN)
    %%% encode only the binary case first
    % morphism generation between two elements or groups of a tensor
    % complex
    %%% note that this is how to map from the image patch 
    % to the midline via the fiber bundle complex.
    % linear regression is used for midline mapping
    % also linear regression is used for SVM via lagrangian and LDA in
    % feature space.
    %%% note that parzen window method is simular to SVM with the addition
    % of eigen vectors computed and then choice made via components in
    % hilbert space rather than summation via the projection 
    % onto the ones-hilbert vector
    %%% note that if the coDomain is a tensor then binary classification 
    % could still be used via the partition of each element of the tensor
    % with a binary relation and a "bit" of wiggle
    %%% note that first attempt will be to have a linear model be the map 
    % and that the linear model will map between features computed from
    % each group of the tensor complex.
    
    
    % step 1: cluster the coDomain and create representation(s)
    % step 
    fidx1 = find(IN.G{IN.cDidx}.IN.D==1);
    for g = 1:size(IN.Didx,1)
        %%%%%%%%%%%%
        % step 1:
        %%%%%%%%%%%%
        % reduce only the sub-cluser of the the gth domain group
        tIN = IN.G{IN.Didx(g)}.IN;
        
        tIN.D = tIN.D(fidx1,:);
        IN.R1{g} = TAC(tIN.D);
        %%%%%%%%%%%%
        % create a hilert rep of the reduction for the subcluster
        tIN.D = {IN.R1{g}.sC};
        tIN.op = 0;
        tIN.nor = 1;
        tIN.P = alpha*ones(1,size(IN.R1{g}.sC,2))/size(IN.R1{g}.sC,2);
        [K tIN] = ker1(tIN);
        [V J1] = eigs(K,IN.G{IN.Didx(g)}.numHil);
        IN.R2{g}.V = V;
        IN.R2{g}.sC = K*V;
        
        %%%%%%%%%%%%
        % step 2:
        %%%%%%%%%%%%
        % reduce the entire gth Domain - via rep type reconstuct
        IN.G{IN.Didx(g)}.op = 3;
        IN.R1{g} = TAC(tIN.D);
        %%%%%%%%%%%%
        % reduce the entire gth Domain - via rep type hilbert
        IN.R2{g}.D{2} = IN.R1{g}.sC;
        %%%%%%%%%%%%
    end
end