function [IN] = stateChange(IN)
%{
%%
% template
%%
IN.D = [ ];
IN.key.op{2}.toOP = [ ]; % numel == ndims(IN.D)
IN.key.op{2}.NC = [ ];  % numel == ndims(IN.D)

%% example _ 1
IN.D = TS;
IN.key.op{2}.toOP = [1 1 1]; % numel == ndims(IN.D);
IN.key.op{2}.NC = [5 5 20];  % numel == ndims(IN.D);
IN.op = 1;
IN = stateChange(IN);
IN.op = 2;
IN = stateChange(IN);

for i = 1:size(TS,3)
    imshow([TS(:,:,i) IN.D(:,:,i)],[]);
    drawnow
end

%% example _ 2
clear IN
IN.D = double(lambda.DOMAIN(:,:,:,lambda.CODOMAIN==1));
IN.key.op{2}.toOP = [1 0 0 0]; % numel == ndims(IN.D);
IN.key.op{2}.NC = [100 5 3 10];  % numel == ndims(IN.D);
% TEST ONE
IN.op = 1;
IN = stateChange(IN);
% TEST TWO
IN.op = 2;
IN = stateChange(IN);
% TEST THREE
% project "new data" set into space
r = round(1 + (size(lambda.DOMAIN,ndims(lambda.DOMAIN))-1).*rand(sum(lambda.CODOMAIN==1),1));
IN.D = double(lambda.DOMAIN(:,:,:,r));
IN.op = 3;
T = stateChange(IN);


%% example _ 3
IN.D = double(squeeze(FA(:,:,1,:,1)));
IN.key.op{2}.toOP = [1 0 0 0]; % numel == ndims(IN.D);
IN.key.op{2}.NC = [100 5 3 10];  % numel == ndims(IN.D);

IN.key.op{2}.dualOP = [0 0 0 0];
IN.op = 1;
IN = stateChange(IN);

IN.key.op{2}.compOP = [];

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
%%%
% 

%}

    if IN.op == 1
        % obtain key set from data and decompose "to" that key set
        IN = of1(IN);
    elseif IN.op == 2
        % from data and key pair - project back to liquid from
        IN = of2(IN);
    elseif IN.op == 3
        IN = of3(IN);
    end
end


% DUAL PROBLEM DONE HERE! - MAY 5 2011
% get key and crystalize
function [IN] = of1(IN)    
    [IN] = operand_op1(IN);
    [IN] = op1(IN,1);
    
    [IN] = operand_op2(IN);
    [IN] = op2(IN,1);    
end

% (re)suspend crystal in solution
function [IN] = of2(IN)    
    [IN] = op2(IN,0);      
    [IN] = op1(IN,0);    
end

% DUAL PROBLEM DONE HERE! - MAY 5 2011
% crystalize given key
function [IN] = of3(IN)        
    [IN] = op1(IN,1);    
    [IN] = op2(IN,1);    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the inverse
function [IN] = operand_op1(IN)
    SZ = size(IN.D);
    tempD = reshape(IN.D,[prod(size(IN.D))/size(IN.D,ndims(IN.D)) size(IN.D,ndims(IN.D))]);
    IN.key.op{1}.operAND.U = mean(tempD,2);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the basis for the crystal
% notes: as i look at this code: i see now that this function is the analog
% of finding the inverse element for the [+] operator. kinda: 
% the inputs of NC make the claim that i "know" something about the inverse
% element of the [*] operator. i am making that claim to move to an
% implementation, however questioning that claim is important.
function [IN] = operand_op2(IN)
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(IN.D)
        SZ = size(IN.D);
        IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))])';
                
        % if the dual should be solved
        if IN.key.op{2}.dualOP(n)
            IN.D = IN.D';
        end
        % if the dual should be solved

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if the key should be obtained!
        if IN.key.op{2}.toOP(n) == 1
            
            [IN.key.op{2}.operAND.U{n} IN.key.op{2}.operAND.BV{n}] = PCA_FIT(IN.D,IN.key.op{2}.NC(n));
            
        else
            
            IN.key.op{2}.operAND.U{n} = zeros(SZ(1));
            IN.key.op{2}.operAND.BV{n} = eye(SZ(1));
            
        end
        % if the key should be obtained!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if the dual should be solved
        if IN.key.op{2}.dualOP(n)
            IN.D = IN.D';
        end
        % if the dual should be solved
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        IN.D = reshape(IN.D,SZ);
        IN.D = permute(IN.D,cvec);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the offset of the crystal
% note the inverse look of sub1 and isub1 - [- > +] and [mean(U) -> U]
function [IN] = op1(IN,direc)        
    SZ = size(IN.D);
    IN.D = reshape(IN.D,[prod(size(IN.D))/size(IN.D,ndims(IN.D)) size(IN.D,ndims(IN.D))]);
    
    for i = 1:size(IN.D,2)
        
        if direc == 0
            IN.D(:,i) = IN.D(:,i) + IN.key.op{1}.operAND.U;
        else
            IN.D(:,i) = IN.D(:,i) - IN.key.op{1}.operAND.U;
        end
        
    end
    
    IN.D = reshape(IN.D,SZ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percipitate into crystal
function [IN] = op2(IN,direc)    
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(IN.D)
        %%% NOTE added wednesday June 29 2011
        if IN.key.op{2}.compOP(n)
        %%% NOTE added wednesday June 29 2011            
        

        

            SZ = size(IN.D);
            IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))]);
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
                simD = (mod1({IN.key.op{2}.operAND.BV{n} newD}));
                IN.key.E{n} = sum((IN.D - simD).^2,1)';
                IN.D = newD;            
            end
            SZ(1) = size(IN.D,1);
            IN.D = reshape(IN.D,SZ);        
        %%% NOTE added wednesday June 29 2011
        end
        %%% NOTE added wednesday June 29 2011
        
        IN.D = permute(IN.D,cvec);
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function call 
function [U BV] = PCA_FIT(M,COM)
    U = mean(M,1);
    for i = 1:size(M,1)
        M(i,:) = M(i,:) - U;
    end    
    COV = cov(M);
    [BV D] = eigs(COV,COM);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decompose function
function [OUT] = decomp(IN)
    

    
    
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
function [K] = genK1(X1)
    K = cov(X1);
end
% KER2 - radial basis funtion
function [K] = genK2(X1,alpha)
    X2 = reshape(X1,[size(X1,1) 1 size(X1,2)]);
    X1 = reshape(X1,[1 size(X1,1) size(X1,2)]);
    K = squeeze(exp(-alpha*sum(bsxfn(@minus,X1,X2).^2,3)));
end
% KER3 - radial basis funtion
function [K] = genK3(X1,X2,alpha)    
    for i = 1:size(X1,1)
        for j = 1:size(X2,1)
            K(i,j) = exp(-alpha*(sum(X1(i,:) - X1(j,:)).^2));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stat_PCA
function [OUT] = statPCA(D,N,R,COM,KERfunc)
    
    cnt = 1;
    for N = 10:10:1000
        for rep = 1:R
            r = round(1 + (size(D,1)-1).*rand(N,1));
            K = KERfunc(D(r,:));

            OUT.U(:,:,rep) = repmat(sum(K,1),[size(K,1) 1]) + repmat(sum(K,2),[1 size(K,1)]) - 2*sum(K(:));
            K = K - OUT.U(:,:,rep);

            %[OUT.BV(:,:,rep),J1] = eigs(K,1);
            %OUT.LAMBA(:,rep) = diag(J1);
            
            [S(cnt,rep,:),J1] = eigs(K,1);
        end
        
        K = KERfunc(squeeze(S(cnt,:,:)));
        DIP = repmat(sum(K,1),[size(K,1) 1]) + repmat(sum(K,2),[1 size(K,1)]) - 2*sum(K(:));
        K = K - DIP;
        [J1,J2] = eigs(K,100);
        VEC(cnt,:) = cumsum(diag(J2)*sum(diag(J2))^-1);
        cnt = cnt + 1;
        N
        %{
        plot(VEC(N-1,:))
        hold on
        
        drawnow
        %}
    end

    
    CLUS = reshape(OUT.BV,[size(OUT.BV,1) size(OUT.BV,2)*size(OUT.BV,3)]);
    CLUS = [CLUS -CLUS];
    
    kidx = kmeans(CLUS',2*COM);
    
    for rep = 1:R
        OUT.BV(:,:,rep)
    end
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRON TWO - pre Kernalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BELOW IS INTRON - WHILE INTERESTING FOR EVOLUTION STUDIES - A TRANSISTION WAS MADE!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%
% 

%}

    if IN.op == 1
        % obtain key set from data and decompose "to" that key set
        IN = of1(IN);
    elseif IN.op == 2
        % from data and key pair - project back to liquid from
        IN = of2(IN);
    elseif IN.op == 3
        
        IN = of3(IN);
    end
end



% get key and crystalize
function [IN] = of1(IN)    
    [IN] = operand_op1(IN);
    [IN] = op1(IN,1);
    
    [IN] = operand_op2(IN);
    [IN] = op2(IN,1);    
end

% suspend
function [IN] = of2(IN)    
    [IN] = op2(IN,0);      
    [IN] = op1(IN,0);    
end

% crystalize given key
function [IN] = of3(IN)        
    [IN] = op1(IN,1);    
    [IN] = op2(IN,1);    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the inverse
function [IN] = operand_op1(IN)
    SZ = size(IN.D);
    tempD = reshape(IN.D,[prod(size(IN.D))/size(IN.D,ndims(IN.D)) size(IN.D,ndims(IN.D))]);
    IN.key.op{1}.operAND.U = mean(tempD,2);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the basis for the crystal
% notes: as i look at this code: i see now that this function is the analog
% of finding the inverse element for the [+] operator. kinda: 
% the inputs of NC make the claim that i "know" something about the inverse
% element of the [*] operator. i am making that claim to move to an
% implementation, however questioning that claim is important.
function [IN] = operand_op2(IN)
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(IN.D)
        SZ = size(IN.D);
        IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))])';
                
        % if the dual should be solved
        if IN.key.op{2}.dualOP(n)
            IN.D = IN.D';
        end
        % if the dual should be solved
        
        % if the key should be obtained!
        if IN.key.op{2}.toOP(n) == 1
            
            [IN.key.op{2}.operAND.U{n} IN.key.op{2}.operAND.BV{n}] = PCA_FIT(IN.D,IN.key.op{2}.NC(n));            
            
            
            
        else
            
            IN.key.op{2}.operAND.U{n} = zeros(SZ(1));
            IN.key.op{2}.operAND.BV{n} = eye(SZ(1));
            
        end
        % if the key should be obtained!
        
        % if the dual should be solved
        if IN.key.op{2}.dualOP(n)
            IN.D = IN.D';
        end
        % if the dual should be solved
        
        IN.D = reshape(IN.D,SZ);
        IN.D = permute(IN.D,cvec);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the offset of the crystal
% note the inverse look of sub1 and isub1 - [- > +] and [mean(U) -> U]
function [IN] = op1(IN,direc)        
    SZ = size(IN.D);
    IN.D = reshape(IN.D,[prod(size(IN.D))/size(IN.D,ndims(IN.D)) size(IN.D,ndims(IN.D))]);
    
    for i = 1:size(IN.D,2)
        if direc == 0
            IN.D(:,i) = IN.D(:,i) + IN.key.op{1}.operAND.U;
        else
            IN.D(:,i) = IN.D(:,i) - IN.key.op{1}.operAND.U;
        end
    end
    
    IN.D = reshape(IN.D,SZ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percipitate into crystal
function [IN] = op2(IN,direc)    
    cvec = 1:ndims(IN.D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(IN.D)
        SZ = size(IN.D);
        IN.D = reshape(IN.D,[SZ(1) prod(SZ(2:end))]);
        % note that here can be inserted a different model AKA the proper
        % mean could be subtracted - this is going be be stripped out into
        % a different function - the REM are for backtracing through the
        % thought process.
        %%%%%%%%%%%%%%%%%%%%%%%% INTRON
        %{
        if direc == 0
            IN.D = IN.key.op{2}.operAND.BV{n}*IN.D;
        else
            IN.D = IN.key.op{2}.operAND.BV{n}'*IN.D;
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%% INTRON
        
        if direc == 0
            % note that "projecting" outwards can be of two types -
            % 1) stochastic or 2) deterministic
            IN.D = mod1({IN.key.op{2}.operAND.BV{n} IN.D});
            %{
            newD = mod1({IN.key.op{2}.operAND.BV{n} IN.D});
            simD = (mod1({IN.key.op{2}.operAND.BV{n}' newD}));
            IN.E{n} = sum((IN.D - simD).^2,1)';
            IN.D = newD;            
            %}
        else
            % generating error upon encryption - error "belongs" to data-key pair            
            newD = mod1({IN.key.op{2}.operAND.BV{n}' IN.D});
            simD = (mod1({IN.key.op{2}.operAND.BV{n} newD}));
            IN.key.E{n} = sum((IN.D - simD).^2,1)';
            IN.D = newD;            
        end
        
        SZ(1) = size(IN.D,1);
        
        IN.D = reshape(IN.D,SZ);        
        IN.D = permute(IN.D,cvec);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function call 
function [U BV] = PCA_FIT(M,COM)
    U = mean(M,1);
    for i = 1:size(M,1)
        M(i,:) = M(i,:) - U;
    end    
    COV = cov(M);
    [BV D] = eigs(COV,COM);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% decompose function
function [OUT] = decomp(IN)
    

    
    
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
% stat_PCA
function [OUT] = statPCA(D,N,R,COM,KERfunc)
    
    for rep = 1:R
        r = round(1 + (size(D,1)-1).*rand(N,1));
        K = KERfunc(D(r,:));
        
        OUT.U(:,:,rep) = repmat(sum(K,1),[size(K,1) 1]) + repmat(sum(K,2),[1 size(K,1)]) - 2*sum(K(:));
        K = K - OUT.U(:,:,rep);
        
        [OUT.BV(:,:,rep),JUNK] = eigs(K,COM);
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bank of kernel functions
function [K] = genK1(X1)
    K = cov(X1);
end

function [K] = genK2(X1,alpha)
    X2 = reshape(X1,[size(X1,1) 1 size(X1,2)]);
    X1 = reshape(X1,[1 size(X1,1) size(X1,2)]);
    K = squeeze(exp(-alpha*sum(bsxfn(@minus,X1,X2).^2,3)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}




%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTRON ONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BELOW IS INTRON - WHILE INTERESTING FOR EVOLUTION STUDIES - A TRANSISTION WAS MADE!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,BV,R] = decomp_ver2()    
    
    [D,key] = sub1(IN.D);
    
    
    BV = sub2(D,NC);
    
    R = sub3(D,BV);
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the offset of the crystal
% note the inverse look of sub1 and isub1 - [- > +] and [mean(U) -> U]
function [D,key] = sub1(D)        
    SZ = size(D);
    D = reshape(D,[prod(size(D))/size(D,ndims(D)) size(D,ndims(D))]);
    
    key.U = mean(D,2);
    
    for i = 1:size(D,2)
        D(:,i) = D(:,i) - key.U;
    end
    
    D = reshape(D,SZ);
end
% and inverse
function [D] = isub1(D,key)    
    SZ = size(D);        
    D = reshape(D,[prod(size(D))/size(D,ndims(D)) size(D,ndims(D))]);

    for i = 1:size(D,2)
        D(:,i) = D(:,i) + key.U;
    end
    
    D = reshape(D,SZ);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the basis for the crystal
% notes: as i look at this code: i see now that this function is the analog
% of finding the inverse element for the [+] operator. kinda: 
% the inputs of NC make the claim that i "know" something about the inverse
% element of the [*] operator. i am making that claim to move to an
% implementation, however questioning that claim is important.
function [key] = sub2(D,NC)
    cvec = 1:ndims(D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(D)
        SZ = size(D);
        D = reshape(D,[SZ(1) prod(SZ(2:end))])';
        
        [key.U{n} key.BV{n}] = PCA_FIT(D,NC(n));
        
        D = reshape(D,SZ);        
        D = permute(D,cvec);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtain the inverse
function [key] = sub1(D)
    SZ = size(D);
    D = reshape(D,[prod(size(D))/size(D,ndims(D)) size(D,ndims(D))]);
    key.U = mean(D,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percipitate into crystal
function [D] = sub3(D,key)    
    cvec = 1:ndims(D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(D)
        SZ = size(D);
        D = reshape(D,[SZ(1) prod(SZ(2:end))]);
        
        D = key.BV{n}'*D;
        
        SZ(1) = size(D,1);
        
        D = reshape(D,SZ);        
        D = permute(D,cvec);
    end
end
% and inverse - suspend
function [D] = isub3(D,key)    
    cvec = 1:ndims(D);
    cvec = circshift(cvec,[0 -1]);
    
    for n = 1:ndims(D)
        SZ = size(D);
        D = reshape(D,[SZ(1) prod(SZ(2:end))]);
        
        D = key.BV{n}*D;
        
        SZ(1) = size(D,1);
        
        D = reshape(D,SZ);        
        D = permute(D,cvec);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ABOVE IS INTRON - WHILE INTERESTING FOR EVOLUTION STUDIES - A TRANSISTION WAS MADE!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}









