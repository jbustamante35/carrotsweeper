function [keep d] = manfredFilter(d)
    d_sub = abs(diff(d,1,2));
    keep = ~any(d_sub > .03,2);
    ratio = mean(d(:,(end-20),:),2).*max(d,[],2).^-1;
    keep = logical(keep .* (ratio > .90));
    d = d(keep,:);
end