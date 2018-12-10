function [m] = initModel(data,num,scale)
    %%%%%%%%%%%%%%    
    % define model space
    MIN = min(data,[],1);
    MAX = max(data,[],1);    
    
    %%%%%%%%%%%%%%
    % make n model - mean
    for n = 1:num
        % iterate over dims
        for d = 1:size(data,2)
            m(n,d) = MIN(d) + (MAX(d)-MIN(d)).*rand(1);
        end
    end
    %%%%%%%%%%%%%%
    % diff scale
    SIGMA = (MAX - MIN)*scale^-1;    
    MIN = -SIGMA;
    MAX = SIGMA;
    %%%%%%%%%%%%%%
    % make n model - var
    for n = 1:num
        % iterate over dims
        for d = 1:size(data,2)
            s(n,d) = MIN(d) + (MAX(d)-MIN(d)).*rand(1);
        end
    end
    m = [m s];
end
