T = rand(10,10,10);
%%%%%
% create vector for folding/unfolding tensor    
cvec = 1:ndims(T);

 for n = 1:ndims(T)        
       
    %%%%%
    % unfolding/folding of tensor
    SZ = size(T);
    T = reshape(T,[SZ(1) prod(SZ(2:end))]);
    uT{n} = mean(T,2);
    
    %T = bsxfun(@minus,T,uT{n});
    T = reshape(T,SZ);   
    
    
    % permute the object
    T = permute(T,cvec);
    
    cvec = circshift(cvec,[0 -1]);
 end
 %%
 %%%%%
% create vector for folding/unfolding tensor    
cvec = 1:ndims(T);


 for n = 1:ndims(T)        
       
    %%%%%
    % unfolding/folding of tensor
    SZ = size(T);
    T = reshape(T,[SZ(1) prod(SZ(2:end))]);
    %uT{n} = mean(T,2);
    
    T = bsxfun(@minus,T,uT{n});
    T = reshape(T,SZ);   
    
    
    % permute the object
    T = permute(T,cvec);
    
    cvec = circshift(cvec,[0 -1]);
 end
 %%
 %%%%%
% create vector for folding/unfolding tensor    
cvec = 1:ndims(T);
cvec = circshift(cvec,[0 -1]);

 for n = 1:ndims(T)        
       
    %%%%%
    % unfolding/folding of tensor
    SZ = size(T);
    T = reshape(T,[SZ(1) prod(SZ(2:end))]);
    uTc{n} = mean(T,2);
    T = reshape(T,SZ);   
    
    
    % permute the object
    T = permute(T,cvec);
 end
%%
%%%%%
 end     
%%%%%
% unfolding/folding of tensor
SZ = size(T);
IN.D = reshape(T,[SZ(1) prod(SZ(2:end))]);