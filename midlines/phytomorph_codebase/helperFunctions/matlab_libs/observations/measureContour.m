function [] = measureClosedContour(C)
    %%% length
    TAN = diff(C,1,1);
    dL = sum(TAN.*TAN,2).^.5;
    L = sum(dL);
    %%% area
    .5*(C(:,2).*TAN(:,1).*

end