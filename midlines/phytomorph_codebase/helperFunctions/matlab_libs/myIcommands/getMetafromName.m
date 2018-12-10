function [] = getMetafromName(N)
    cmd = ['imeta ls -d ' N];
    [status query] = system(cmd);
    query
end