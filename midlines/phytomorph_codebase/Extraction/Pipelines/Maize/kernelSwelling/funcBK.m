function [c] = funcBK(a,b,x)
    %c = a*(1-exp(-b*x));
    c = bsxfun(@times,a,(1-exp(-b*x)));
    %c = b*x.*(a + x).^-1;
end

%%
%{
t = 0:150;
kvec = linspace(.01,.02,10);
max = .11;
for k = 1:numel(kvec)
    d = func(max,kvec(k),t);
    [J midx] = min((abs(d) - max/2));
    dd = gradient(d);
    n(k) = dd(midx);
    plot(d)
    hold all
end
figure;
plot(kvec,n);
%}