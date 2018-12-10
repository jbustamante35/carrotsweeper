function [S] = getPCAdy(X,IDX)
    for tm = 1:size(X,2)
        [S(tm).S S(tm).C S(tm).U S(tm).E S(tm).L S(tm).ERR S(tm).LAM] = PCA_FIT_FULL([X(:,IDX),X(:,tm)],3);
        fidx = find(S(tm).LAM < eps);
        S(tm).LAM(fidx) = 0;
        if tm > 1
            for d = 1:size(S(tm).E,2)
                if S(tm-1).E(:,d)'*S(tm).E(:,d) < 0
                    S(tm).E(:,d) = -S(tm).E(:,d); 
                end
            end
        end
    end
end