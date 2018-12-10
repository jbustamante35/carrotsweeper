function [P] = deathCycle(FNp,P,percent)
    if numel(P) ~= 0
        % get number of ind to check
        sidx = randi(numel(P),round(percent*numel(P)),1);
        
        if numel(P)*percent < .5
            sidx = randi(numel(P),1);
        end
        
        grp = P(sidx);
        fprintf(['Selected-death:' num2str(numel(grp)) ':' num2str(numel(P)) '\n'])
        % loop over each ind
        parfor e = 1:numel(grp)
            dead = logical(0);
            % get indivdual
            ind = grp(e);
            % get chR number
            UQ = unique(ind.cn(:,1));
            % loop over each ch
            for u = 1:numel(UQ)
                % get ch index
                cidx = ind.cn(:,1) == UQ(u);
                % get chromosome
                ch = ind.ch(cidx);
                % check death sites
                dead = dead | any(round(ch(FNp(u).sites)) == FNp(u).values);
                if dead;break;end
            end
            rm(e) = logical(dead);
        end
        P(rm) = [];
    end
end