function [path] = ipwd()
    [status, path] = system('ipwd');
    path = path(1:end-1);
end