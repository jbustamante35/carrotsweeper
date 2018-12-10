function [rf] = recombinationFrequency(md)
    rf = ones(1,size(md,1));
    rf = diag(rf);
    for i = 1:size(md,1)
        for j = (i+1):size(md,1)
            fprintf(['starting ibd:' num2str(i) ':' num2str(j) '\n']);
            tic
            rf(i,j) = calc(md(i,:),md(j,:));
            %rf(j,i) = rf(i,j);
            fprintf(['ending ibd:' num2str(i) ':' num2str(j) ':' num2str(toc) '\n']);
        end
    end
    rf = rf + triu(rf)';
end

function [rf] = calc(v1,v2)
    UQ1 = unique(v1);
    UQ2 = unique(v2);
    for u1 = 1:numel(UQ1)
        tmp1 = v1 == UQ1(u1);
        for u2 = 1:numel(UQ2)
            tmp2 = v2 == UQ2(u2);
            rf(u1,u2) = sum(tmp1.*tmp2)*numel(tmp1)^-1;
        end
    end
    rf;
end