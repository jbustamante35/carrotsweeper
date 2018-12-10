function [P] = birthCycle(CHp,P,percent,xM1,xM2)
    try
        % get number of ind to mate
        sidx1 = randi(numel(P),round(percent*numel(P)),1);
        sidx2 = randi(numel(P),round(percent*numel(P)),1);

        if numel(P)*percent < .5
            sidx1 = randi(numel(P),1);
            sidx2 = randi(numel(P),1);
        end

        grp1 = P(sidx1);
        grp2 = P(sidx2);
        
        
        fprintf(['Selected-birth:' num2str(numel(grp1)) ':' num2str(numel(P)) '\n']);
        parfor e = 1:numel(grp1)
            if nargin == 3
                [Pn(e)] = crossInv(grp1(e),grp2(e),CHp);
            else
                [Pn(e)] = crossInv(grp1(e),grp2(e),CHp,xM1,xM2);
            end
        end
        P = [P Pn];
    catch ME
        ME
    end
end