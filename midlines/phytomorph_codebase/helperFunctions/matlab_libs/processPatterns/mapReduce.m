function [] = mapReduce(F,T,X)
    % select out the nodes from store X
    fidx = find(X.queryNodeProp(T));
    % operate on each element of X
    for e = 1:numel(fidx)                
        % select current node
        x = X.d{fidx(e)};
        % apply F to x (element of) X
        % F also will store and link - bound from above
        F(x);
        % clear conNode
        clear x
    end
end