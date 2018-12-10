function [y] = myOp0(m,v,k)
    y = mtimesx(v(k,:),m);
    %y = v(k,:)*m;
end