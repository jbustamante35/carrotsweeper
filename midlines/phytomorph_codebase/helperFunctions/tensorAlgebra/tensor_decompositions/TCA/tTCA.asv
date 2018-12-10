function [IN] = tTCA(IN)
%{
%%%%%%%%%%%%%%%%
% template
%%%%%%%%%%%%%%%%

IN.D = []; data reduction - last
IN.key.op{2}.toOP = [1 0 0 0]; % numel == ndims(IN.D);
IN.key.op{2}.NC = [100 5 3 10];  % numel == ndims(IN.D);
IN.key.op{2}.compOP = [];
IN.op = 1;
IN = stateChange(IN);

%}

%{
%%%%%%%%
% notes: due to the observation of symmetry in the code: [+,*] are being defined below:
% the operators [+,*] are admitting/displaying both mono and bi valent
% behavior. when mono valent - the second element is derived from the
% first, as in "the mean" or expected value operation (in the case of [+]
% operator) AND "eigen" style operator (in the case of [*]).
%%%
% in some way-we have a group/ring style structure emerging. that is we
% have two operators and a "set" being the vectors/tensors of data.  
%%%
% in addition to the algebric struction: i am looking for the analogy for the basis vectors/tensors
% these seem to place footing to the set and help to form a ring or a
% field.  by this i mean that there may not be an inverse "function"
%%%
% i am also in the search of the "key" in the below code.  it serves, at
% some level, to help by 1) seperating a "new" data trial into precipitate
% and supernate 2) providing structure to the idea of ring vs field.  that
% is: while a "new trial" may not be in the precipitate and there for might
% belong to a ring like structure - ONLY if projected to the basis - it
% also can be described as a field if the data is not LOST.
%%%%%%%%
% notes: added 4/20 : reconstruction error added - there are two classes of
% methodologies for determining the DIM of the reduction - rule based OR
% defined - rule based reduction will give rise to the defined whereas the
% defined type will start with DIM priori
%%%%%%%%
% notes: added 11.12.07: bookend the tensor decomposition with the
% linearized version.  unwrap the compressed tensor and decompose it too.
% added a flag to have this or not.

% 
%}

    if IN.op == 1
        % obtain key set from data and decompose "to" that key set
        % return the decompostion and the key that generate it
        % note the sequece :
        % data -> representation <-> reconstruction + err        
        IN = of1(IN);
    elseif IN.op == 2
        % from data and key pair - project back to liquid from
        % note that there are a few ways to look at this.
        % encryption, crystals, and language (high level abstraction(s))
        IN = of2(IN);
    elseif IN.op == 3
        % from data and key pair - encode/encrypt,crystalize,reduce
        IN = of3(IN);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get key and crystalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IN] = of1(IN)    
    [IN] = operand_op0(IN);
    [IN] = op0(IN,1);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (re)suspend crystal in solution
function [IN] = of2(IN)    
    [IN] = op2(IN,0);      
    [IN] = op1(IN,0);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crystalize given key
function [IN] = of3(IN)        
    [IN] = op1(IN,1); 
    [IN] = op2(IN,1);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the basis for the crystal
% notes: as i look at this code: i see now that this function is the analog
% of finding the inverse element for the [+] operator. kinda: 
% the inputs of NC make the claim that i "know" something about the inverse
% element of the [*] operator. i am making that claim to move to an
% implementation, however questioning that claim is important.
% NOTES (added.11.12.06): keeping the assumption about the inverse of [+],
% but recalling that magma, ring, field, have something to do with the
% future of this.
function [IN] = operand_op0(IN)
    %%%%%
    % create vector for folding/unfolding tensor
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(IN.D)
        %%%%%
        % unfolding/folding of tensor
        SZ = size(IN.D);
        IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))])';
        
        %%%%%
        % if operation should be performed
        %if IN.key.op{2}.toOP(n) == 1
        if IN.key.op{2}.NC(n) > 0            
            %%%%%
            % return the operand and basis vectors
            [IN.key.op{2}.operAND.U{n} IN.key.op{2}.operAND.BV{n}] = PCA_FIT(IN.D,IN.key.op{2}.NC(n));
        else
            %%%%%    
            % return "nothing" elements
            % return the iD elements of the algebra
            IN.key.op{2}.operAND.U{n} = 0;
            IN.key.op{2}.operAND.BV{n} = 1;
        end
        
        %%%%%
        % inverse fold
        IN.D = reshape(IN.D,SZ);
        IN.D = permute(IN.D,cvec);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percipitate into crystal
% use the operaANDS on the objects
function [IN] = op2(IN,direc)    
    %%%%%
    % create vector for folding/unfolding tensor    
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    %%%%%
    % for each "dim" 
    for n = 1:ndims(IN.D)        
        %%%%%
        % if a computation on the obect is being requested
        if IN.key.op{2}.NC(n) > 0
            %%%%%
            % unfolding/folding of tensor
            SZ = size(IN.D);
            IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))]);
            
            %%%%%
            % use the model for the data
            % note that here can be inserted a different model AKA the proper
            % mean could be subtracted - this is going be be stripped out into
            % a different function - the REM are for backtracing through the
            % thought process.
            if direc == 0
                % note that "projecting" outwards can be of two types -
                % 1) stochastic or 2) deterministic
                IN.D = mod1({IN.key.op{2}.operAND.BV{n} IN.D});            
            else
                % generating error upon encryption - error "belongs" to data-key pair            
                newD = mod1({IN.key.op{2}.operAND.BV{n}' IN.D});
                IN.D = newD;            
            end
            
            %%%%%
            % inverse fold
            SZ(1) = size(IN.D,1);
            IN.D = reshape(IN.D,SZ);                
        end
        %%%%%
        % permute the object
        IN.D = permute(IN.D,cvec);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function call 
function [U BV] = mod_f0(M,c,n)
    
    U = cop(@mean,M,n);    
    M = M - U;
    
    [BV D] = eigs(COV,COM);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "model function"-- without mean considered
function [OUT] = mod1(IN)
    OUT = IN{1}*IN{2};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "model function"-- with mean considered
function [OUT] = mod2(IN)
    OUT = IN{1}*IN{2} + IN{3};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bank of kernel functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KER1 - inner product - COV
function [IN] = ker(IN)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % gererate kernel
    switch IN.type
        case 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % covar kernel
            K = cov(IN.D{1});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% select which D from the cell
            idx = size(IN.D,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% allocate the kernel
            K = zeros(size(IN.D{idx},1),size(IN.D{1},1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% calc the kernel 
            for i = 1:size(IN.D{idx},1)
                for j = 1:size(IN.D{1},1)
                    K(i,j) = mvnpdf(IN.D{1}(j,:),IN.D{idx}(i,:),IN.P);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% normalize the kernel
    if IN.nor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% obtain offsets
        cU = mean(K,2);
        if IN.op == 1
            rU = IN.rU;
            oU = IN.oU;
        else
            rU = mean(K,2);
            oU = mean(K(:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% store
        if IN.op == 0
            IN.rU = rU;
            IN.oU = oU;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% calc and normalize - this is a bit slower
        %%% speed up later
        for i = 1:size(IN.D{idx},1)
            for j = 1:size(IN.D{1},1)
                K(i,j) = K(i,j) - rU(j) - cU(i) + oU;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% return kernel
    IN.K = K;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

