function [ibdM] = idb(C1,C2)
    ibdM = 0;
    UQ = unique(C1.cn(:,1));
    for u = 1:numel(UQ)
        cidx = C1.cn(:,1) == UQ(u);
        ch1 = C1.ch(cidx);
        ch2 = C2.ch(cidx);
        ch1 = foldCH(ch1);
        ch2 = foldCH(ch2);
        CHt = [ch1,ch2];
        parfor k = 1:size(CHt,1)
            % find number of potnetial alleles
            tmp = unique(CHt(k,:));
            tmpN = [];
            % find number of occurances of each allel - allel freq?
            for u2 = 1:numel(tmp)
                tmpN(u2) = sum(CHt(k,:) == tmp(u2));
            end
            % find non identcal
            ibdM(k) = sum(tmpN ~= 1)/4;
        end
    end
    ibdM = mean(ibdM);
end