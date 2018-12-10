function [curveBank] = selectClosedCurves(curveBank)
    rm = [];
    for e = 1:numel(curveBank)
        if ~all(curveBank(e).data(:,1) == curveBank(e).data(:,end))
            rm(e) = 1;
        end
    end
    curveBank(find(rm)) = [];
end