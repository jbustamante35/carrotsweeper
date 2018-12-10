function [p] = LOGprobData(data,model)
    % non-rotational
    % model [u,s]
    dsz = size(data,2);
    idx = (dsz+1);
    % get mu and var from model
    mu = model(1:dsz);
    var = model(idx:idx+(dsz-1));
    % eval the log
    delta = bsxfun(@minus,data,mu);
    delta = bsxfun(@times,delta,var.^-2);
    p = sum(delta.*delta,2).^.5;
end