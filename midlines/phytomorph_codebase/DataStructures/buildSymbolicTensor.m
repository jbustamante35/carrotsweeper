function [F] = buildSymbolicTensor(orderVec)
    syms F
    T = sym(['T_'], orderVec,'real');
    t = sym(1);
    newargs = [];
    for e = 1:numel(orderVec)
        V = sym(['V' num2str(e) '_%d'],[orderVec(e) 1],'real');
        t = (t*V');
        t = t(:);
        newargs = [newargs;V];
    end
    lambda = sym('l',[numel(t) 1]);
    F(newargs) = T(:)'*t(:);
    %F(newargs) = subs(F,lambda,t)
end

%{  
    [T] = buildSymbolicTensor([4 5 6]);
    

%}