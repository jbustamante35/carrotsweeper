function [sM] = cosegregationMetric(md)
    sM = .5*ones(1,size(md,1));
    sM = diag(sM);
    for i = 1:size(md,1)
        tmpV = md(i,:);
        parfor j = (i+1):size(md,1)
            fprintf(['starting csm:' num2str(i) ':' num2str(j) '\n']);
            tic
            csd = calc(tmpV,md(j,:));
            sM(i,j) = nanmean(log(csd(:)).*csd(:));
            %sM(i,j) = calc(md(i,:),md(j,:));
            %sM(j,i) = sM(i,j);
            fprintf(['ending csm:' num2str(i) ':' num2str(j) ':' num2str(toc) '\n']);
        end
    end
    sM = sM + triu(sM)';
end

function [csd] = calc(v1,v2)
    % extract all states from data - should declare all states
    UQ1 = unique(v1);
    UQ2 = unique(v2);
    for u1 = 1:numel(UQ1)
        tmp1 = v1 == UQ1(u1);
        for u2 = 1:numel(UQ2)
            tmp2 = v2 == UQ2(u2);
            csd(u1,u2) = sum(tmp1.*tmp2)*numel(tmp1)^-1;
        end
    end
end