function [mapper] = masterMap(X,Y,para,type)
    
    if nargin == 1
        % return default para for 
    end
    
    % init default parameters
    switch type
        case 'plsregression'
            para.MaxComp = min(size(X,1)-1,size(X,2));
            para.MinComp = 1;
    end
    
    switch type
        case 'plsregression'
            cnt = 1;
            for n = para.MinComp:para.MaxComp
                trainFunc = @(X,Y)myPLSwrapper(X,Y,n);
                vals{cnt} = crossval(trainFunc,X,Y);
            end
        case 'cca'
        case 'cnn'
        case 'pcr'
        case 'ols'
    end
    
end

function [metric] = myPLSwrapper(X,Y,n)
    [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(X,Y,n);
    metric = MSE(2,end);
end



%{




%}