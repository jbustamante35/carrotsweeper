function [t] = multv(t,v)
    % each v is [1 N] and s = [N s];
    %t.s = [numel(v) t.s];
    %t.d = kron(t.d,v);
    t.s = [t.s numel(v)];
    t.d = kron(v,t.d);
end

%{
    v1 = [0 0 1];
    v2 = [1 0];
    v3 = [1 1];
    t.d = v1;
    t.s = numel(v1);
    t = multv(t,v2);
    a = reshape(t.d,t.s);
    t = multv(t,v3);
    b = reshape(t.d,t.s);
%}