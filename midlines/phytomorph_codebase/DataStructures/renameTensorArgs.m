function [nT] = renameTensorArgs(T,orderVec,newBase)
    syms nT;
    newargs = [];
    for order = 1:numel(orderVec)
        V = sym(['V' num2str(order) '_%d'],[orderVec(order) 1]);
        nB = sym([newBase(order) '_%d'],[orderVec(order) 1]);
        for e = 1:numel(V)
            T = subs(T,V(e),nB(e));
        end
        newargs = [newargs;nB];
    end
    nT(newargs) = T;
end

%{
    syms alpha theta rho
    tC_rot = [cos(alpha);sin(alpha);-sin(alpha);cos(alpha)];
    
    tC_polar2car = [cos(theta);0;0;sin(theta)];

    [T] = buildSymbolicTensor([2 2]);
    

    rot = fillTensorComponents(T,tC_rot,[2 2]);
    rot = renameTensorArgs(rot,[2 2],['O' 'I']);


    p2c = fillTensorComponents(T,tC_polar2car,[2 2]);
    p2c = renameTensorArgs(p2c,[2 2],['O' 'I']);


    syms x y I_1 I_2 O_1 O_2
    %rot(1,0,x,y)

    g(O_1,O_2,I_1,I_2) = compose(rot,p2c,[I_1],[O_1]);
    g(O_1,O_2,I_1,I_2) = compose(g,p2c,[I_2],[O_2]);
    g(1,0,rho,rho)
    G = [g(1,0,I_1,I_2);g(0,1,I_1,I_2)];











%}