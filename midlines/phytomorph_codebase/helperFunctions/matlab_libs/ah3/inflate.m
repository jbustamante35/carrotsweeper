function [model] = inflate(model,d,target,groups)
    dist_min = ah3(model,d,target,groups);

    dsz = size(d,2);
    idx = (dsz+1);
    mu = model(1:dsz);
    var = model(idx:idx+(dsz-1));
    max_expand_percent = 10;
    count_expand_threshold = .1;
    expandV = linspace(1,max_expand_percent,100);
    
    for se = 1:numel(var)
        for e = 1:numel(expandV)
            n_var = var;
            n_var(se) = expandV(e)*n_var(se);
            n_model = [mu n_var];
            dist(se,e) = ah3(n_model,d,target,groups);
        end
    end
    
    select_for_expand = find(dist(:,end) < dist_min*(count_expand_threshold + 1));
    
    n_var = var;
    for i = 1:numel(select_for_expand)
        n_var(select_for_expand(i)) = max_expand_percent*n_var(select_for_expand(i));
    end
    
    
    n_model = [mu n_var];
end