function [t] = nmmult(t,m,n)
    % create dimension vector - swtich order of 1,n
    dm = 1:ndims(t.d);
    dm(1) = n;
    dm(n) = 1;
    % premute to focus on dim(n)
    t.d = permute(t.d,dm);
    % reshape for multiplication and obtain size
    [t.d s] = r(t.d);
    % perform multiplication
    t.d = t.d*m';
    % refresh size
    s(1) = size(t.d,2); % alternate s(1) = size(m,1);
    % inverse reshape
    t.d = ir(t.d,s);
    % inverse permute
    t.d = ipermute(t.d,dm);
end



%%%%%%%%%%%%%%%%
% reshape
function [d s] = r(d)
    s = size(d);    
    d = reshape(d,[s(1) prod(s(2:end))])';
end

%%%%%%%%%%%%%%%%
% inverse reshape
function [d s] = ir(d,s)    
    d = reshape(d',s);
end