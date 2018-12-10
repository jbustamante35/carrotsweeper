function [S] = sampleG(I,N,G)
    S = zeros(size(G,3),size(N,2));
    for i = 1:size(G,3)
        No = G(:,:,i)*N;
        S(i,:) = ba_interp2(I,No(2,:),No(1,:));
    end    
end