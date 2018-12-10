function [wFunc] = getPNN_func(X,Y,d)
    tmp = tempname;
    tmp(1:5) = [];
    functionName = ['/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/nets/' tmp '.m'];
    [p,nm,ext] = fileparts(functionName);
    
    
    [~,nX,U,E] = PCA_FIT_FULL_T(X,20);
    
    
    for e = 1:10
        net{e} = patternnet(d);
        %Y = [Y==0;Y==1];
        [net{e},tr{e}] = train(net{e},nX,Y,'useParallel','yes');
    end
    clear Y clear X nX
    for e = 1:10
        bb(e) = tr{e}.best_perf;
    end
    [~,midx] = min(bb);
    %funcObject.net = net{midx};
    
    %{
    %%%%%
    G = net{midx}(nX);;
    Y = zeros(size(Y));
    Y(G>.5) = 1;
    
    
    for e = 1:10
        net{e} = patternnet(d);
        %Y = [Y==0;Y==1];
        [net{e},tr{e}] = train(net{e},nX,Y,'useParallel','yes');
    end
    
    for e = 1:10
        bb(e) = tr{e}.best_perf;
    end
    [~,midx] = min(bb);
    funcObject.net = net{midx};
    %%%%
    %}
    
    genFunction(net{midx},functionName);
    clear net
    funcObject.func = str2func(nm);
    
    
    wFunc = @(X,e0,e1)funcObject.func(PCA_REPROJ_T(X,E,U));
    
    
    %pF = partialFunction(func,'maizeSeedlings');
    %pF.publish();
    %{
    Weights = mynLDA(X',Y',1,d);
    U = mean(X,2);
    wFunc = @(X,e0,e1)(((Weights')*bsxfun(@minus,X,U)));
    %}
end