function [T] = fillTensorComponents(T,values,orderVec)
    tmp = sym(['T_'], orderVec);
    tmp = tmp(:);
    for e = 1:numel(tmp)
        T = subs(T,tmp(e),values(e));
    end
end

%{
    syms alpha theta rho
    tC_rot = [cos(alpha);sin(alpha);-sin(alpha);cos(alpha)];
    
    tC_polar2car = [cos(theta);0;0;sin(theta)];

    [T] = buildSymbolicTensor([2 2]);
    
    rot = fillTensorComponents(T,tC_rot,[2 2]);
    p2c = fillTensorComponents(T,tC_polar2car,[2 2]);
    syms x y
    rot(1,0,x,y)

    compose(rot,p2c,


%}