function [paraPRE_X,paraNPRE_X,paraPRE_Y,preData,postData] = extractTraits(D,N,preX_model,norX_model,preY_model)

    preData = D(1:N,:);
    postData = D((N+1):end,:);
    
    paraPRE_X = preX_model.func(preData);
    
    paraNPRE_X = norX_model.func(preData);
    
    paraPRE_Y = preY_model.func(preData,postData);
    
end