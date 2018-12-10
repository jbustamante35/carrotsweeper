function [ld] = calcLD(P,type)
    CH = zeros(size(P(1).ch,1),numel(P));
    for e = 1:numel(P)
        CH(:,e) = P(e).ch;
    end
    CH = round(CH);
    for e1 = 1:size(CH,1)
        v1 = CH(e1,:);
        parfor e2 = (e1+1):size(CH,1)
            fprintf(['starting ld:' num2str(e1) ':' num2str(e2) '\n']);
            tic
            ld = calc(v1,CH(e2,:));
            fprintf(['ending ld:' num2str(e1) ':' num2str(e2) ':' num2str(toc) '\n']);
        end
    end
end

function [ld] = calc(v1,v2)
    % lower level function
    % v1 := marker reads at loci1
    % v2 := marker reads at loci2
    uq1 = unique(v1);
    uq2 = unique(v2);
    for u = 1:numel(uq1)
        p1(u) = sum(v1==uq1(u))/numel(v1);
    end
    for u = 1:numel(uq2)
        p2(u) = sum(v2==uq2(u))/numel(v2);
    end
    
    pti = p1'*p2;
    
    for u1 = 1:numel(uq1)
        tmp1 = v1 == uq1(u1);
        for u2 = 1:numel(uq2)
            tmp2 = v2 == uq2(u2);
            tmp = sum(tmp1 & tmp2);
            pt(u1,u2) = tmp/numel(v1);
        end
    end
    
    ld = pt - pti;
end