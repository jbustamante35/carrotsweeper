classdef qO < handle
    properties
        func;
    end
    
    
    methods
        function [obj] = qO(func)
            obj.func{1}{1} = func;
        end
        
        function [r] = mtimes(obj,v)
            r = sv(size(v,1),size(v,2),numel(obj.func));
            % if v is sparse
            if issparse(v)
                % for each map
                for m = 1:numel(obj.func)
                    % for each column vector
                    for e = 1:size(v,2)
                        % find the support domain
                        fidx = find(v(:,e));
                        % create sparse range
                        tmp = sparse(size(v,1),1);
                        % for each function within a map - supports
                        % operator superposition
                        for f = 1:numel(obj.func{m})
                            % get the range support
                            nidx = obj.func{m}{f}{1}(fidx);
                            values = obj.func{m}{f}{2}(fidx);
                            tmp(nidx) = tmp(nidx) + v(fidx,e).*values;
                        end
                        ridx = find(tmp);
                        r(ridx,e,m) = tmp(nidx);
                    end
                end
            end
        end
    end
end

%{
    perTime = 1;
    q1 = qO(@(C)(C+perTime));



    % position operator
    P = qO({@(X)(X),@(X)(X)});
    % position operator test
    p = sv(1000,1);
    p(10,1) = .5.^.5;
    p(12,1) = .5.^.5;
    q = P*p;
    e = p*(P*p);



    v = sv(10^9,1);
    %v = sparse();
    v(1,1) = 1;
    tic
    w = q1*v;
    toc

    prob = qO(@(C)(C*fraction));
    




%}