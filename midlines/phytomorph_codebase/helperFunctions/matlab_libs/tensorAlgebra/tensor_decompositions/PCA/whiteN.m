function [M,U,S] = whiteN(M)
    U = mean(M,1);
    for i = 1:size(M,1)
        M(i,:) = M(i,:) - U;
    end
    S = std(M,1,1);
    for i = 1:size(M,1)
        M(i,:) = M(i,:).*S.^-1;
    end
end