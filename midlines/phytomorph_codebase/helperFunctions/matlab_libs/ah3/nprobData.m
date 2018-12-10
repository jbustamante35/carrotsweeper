function [p] = nprobData(data,model)
    %%%%%%%%%%%%%%%%%%%
    % noramlized prob of data
    %%%%%%%%%%%%%%%%%%%
    % for each data
    dm = size(data,2);
    % get prob of data
    p = probData(data,model);
    % get max of model
    mx = probData(model(1:dm),model);
    % normalize to one
    p = p/mx;
end