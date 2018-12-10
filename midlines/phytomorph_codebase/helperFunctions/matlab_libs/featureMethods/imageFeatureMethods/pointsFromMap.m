function [pointList] = pointsFromMap(featureMap,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run non-max suppresion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rad = para.nonmaxsuppts.rad.value;
    lmM = nonmaxsuppts(featureMap, rad);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold = para.binaryOperator.threshold.value;
    operator = para.binaryOperator.op.value;
    bcM = operator(featureMap,threshold);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % logical call for feature point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logicMap = cat(3,bcM,lmM);
    logicMap = attachMask(logicMap,para.MSK);
    pointList = featureSelection(logicMap);
end


% this needs to be better
% if the mask is not present then do not attach
% if the mask does not fit then do not attach
% if the mask fits then attach
function [logicMap] = attachMask(logicMap,MSK)
    if ~isempty(MSK)
        try
            logicMap = cat(3,logicMap,MSK);
        catch
            return;
        end
    end
end