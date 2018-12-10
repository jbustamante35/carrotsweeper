function [target] = expectedKernelNumber(fileName)
    try
    % find parts of file(s)
    [p n e] = fileparts(fileName);
    % remove the name
    n(end-4:end) = [];
    % find the key char
    fidx = strfind(n,'_');
    %
    n(1:fidx(2)) = [];
    n = ['_' n '_'];
    fidx = strfind(n,'_');
    target = numel(fidx) - 1;
    catch ME
        ME
        target = 0;
    end
end