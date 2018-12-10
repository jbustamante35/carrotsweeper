function [s] = genString(A,L,D)
    s = '';
    if L > 0
        r = randi(numel(A),1,L);
        for e = 1:L
            s = [s,A{r(e)},D];
        end
        s(end) = [];
    end
end

%{
    A = {'1','2','3','4','5','6','7','8','9'};
    D = ',';
    L = 5;
    s = genString(A,L,D);
%}