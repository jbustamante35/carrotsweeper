function [ A ] = skel2adj( S )
%SKEL2ADJ creates a sparse adjacency matrix to be used for Djikstra's alg.
%   
S = padarray(S, [1 1], 0);

[r,c] = find(S);
A = zeros(size(r,1));
for i = 1:size(r,1)
    Y = r(i); X = c(i);
    t = S((Y-1):(Y+1),(X-1):(X+1)); % Should padarray(S) to prevent errors at edges
    t(2,2) = 0;
    [y,x] = find(t);
    y = y - 2; x = x - 2;
    
    connect = zeros(size(y));
    for j=1:size(connect,1);
        connect(j,1) = find( r == (Y + y(j,1)) & c == (X + x(j,1)));
    end
    
    A(i, connect) = 1;
    %i
end
    
A = sparse(A);

end

