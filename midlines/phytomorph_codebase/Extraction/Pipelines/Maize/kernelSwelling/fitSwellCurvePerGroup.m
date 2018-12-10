function [fixed,random,stats] = fitSwellCurvePerGroup(D)
    % time for curves
    TM = 1:size(D,2);
    TM = TM - 1;
    TM = repmat(TM,[size(D,1) 1]);
    groups = (1:size(D,1))';
    groups = repmat(groups,[1 size(D,2)]);
    options = statset('Display','iter');
    ridx = isnan(D(:));
    TM(ridx) = [];
    groups(ridx) = [];
    D(ridx) = [];
    finalGuess = mean(D(:,end));
    [fixed,PSI1,stats,random] = nlmefit(TM(:),D(:),groups(:),[],@(alpha,t)swellCurveFit(alpha,t),[finalGuess .01],'Options',options);  
end