function [C] = fuseGam(G0,G1)
    UQ = unique(G0.cn(:,1));
    C.cn = [];
    C.ch = [];
    for u = 1:numel(UQ)
        p0 = G0.cn(:,1) == UQ(u);
        p1 = G1.cn(:,1) == UQ(u);
        
        C.ch = [C.ch;[G0.ch(p0);G1.ch(p1)]];
        C.cn = [C.cn;[UQ(u)*ones(sum(p0) + sum(p1),1),[zeros(sum(p0),1);ones(sum(p1),1)]]];
    end
end