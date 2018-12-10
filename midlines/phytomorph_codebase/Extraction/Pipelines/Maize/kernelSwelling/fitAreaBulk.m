function [p ep er] = fitAreaBulk(d,ft)
    for tr = 1:size(d,2)        
        %[p(tr,:) er(tr)] = fminsearch(@(X)mySwellFit(d(:,tr),X),[10^4 .1]);
        %[p(tr,:) er(tr)] = fminunc(@(X)mySwellFit(d(:,tr),X),[10^4 .1]);
        
        [p(tr,:) er(tr)] = fminsearch(@(X)mySwellFit(d(:,tr)',X),[.2 .001]);
        %[p(tr,:) er(tr)] = fminunc(@(X)mySwellFit(d(:,tr)',X),[1 1]);
        ep(tr,:) = func(p(tr,1),p(tr,2),0:(ft-1));
        tr
    end
end