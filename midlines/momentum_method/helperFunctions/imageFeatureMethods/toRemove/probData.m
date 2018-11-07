function [p] = probData(data,model)
    % non-rotational
    % model [u,s]
    d = size(data,2);
    idx = (d+1);
    mu = model(1:d);
    var = model(idx:idx+(d-1));
    p = mvnpdf(data,mu,abs(var).^2);
end



