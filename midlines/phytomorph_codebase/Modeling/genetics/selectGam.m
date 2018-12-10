function [G0,G1] = selectGam(C)
    G0.ch = [];
    G0.cn = [];
    G1.ch = [];
    G1.cn = [];
    UQ = unique(C.cn(:,1));
    for u = 1:numel(UQ)
        subUQ = find(C.cn(:,1) == UQ(u));
        p0 = find(C.cn(subUQ,2) == 0);
        p1 = find(C.cn(subUQ,2) == 1);
        G0.ch = [G0.ch;C.ch(subUQ(p0))];
        G1.ch = [G1.ch;C.ch(subUQ(p1))];
        G0.cn = [G0.cn;UQ(u)*ones(numel(p0),1)];
        G1.cn = [G1.cn;UQ(u)*ones(numel(p1),1)];
    end
end