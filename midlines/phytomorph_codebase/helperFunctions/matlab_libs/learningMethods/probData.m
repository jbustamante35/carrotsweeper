function [p] = probData(data,model)
    
    % non-rotational
    % model [u,s]
    d = size(data,2);
    idx = (d+1);
    mu = model(1:d);
    var = model(idx:idx+(d-1));
    
    p = mvnpdf(data,mu,var);
    
    %alpha = 100;
    %data = data*alpha;
    %mu = mu*alpha;
    %var = var*alpha;
    %p = mvnpdf(data,mu,var.^2);
    
    
    
    %{
    alpha = 1000;
    p = bsxfun(@minus, data,mu);
    p = bsxfun(@times,p,var.^-1);
    p = sum(p.*p,2).^.5;
    p = p * alpha^-1;
    p = normpdf(p);
    %}
end



