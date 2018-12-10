function [k] = kinship(P)
    k = ones(1,numel(P));
    k = diag(k);
    for i = 1:numel(P)
        parfor j = (i+1):numel(P)
            fprintf(['starting ibd:' num2str(i) ':' num2str(j) '\n']);
            tic
            k(i,j) = ibd(P(i),P(j));
            %k(j,i) = k(i,j);
            fprintf(['ending ibd:' num2str(i) ':' num2str(j) ':' num2str(toc) '\n']);
        end
    end
    k = k + triu(k)';
end
