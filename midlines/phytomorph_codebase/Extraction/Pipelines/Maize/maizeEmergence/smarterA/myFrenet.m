function [threadMetric] = myFrenet(sequence,lambdaVec,funcG)
    % try a windowed approach
    threadMetric = [];
    for scale = 1:numel(lambdaVec)
        lambda = lambdaVec(scale);
    
        % half window size
        h = (lambda-1)/2;
        % pad array by amount
        tmpPADsequence = padarray(sequence,[0 h],'replicate','both');
        Nsequence = [];
        % apply func to each window
        for e = 1:(size(tmpPADsequence,2)-(2*h))
            Nsequence(:,e) = funcG(tmpPADsequence(:,e:(e+2*h)));
        end

        % flag for upper vs diag
        upper = 0;
        % make lambda the number of dims of thread - curve
        lambda = size(Nsequence,1);
        % make half window size
        h = (lambda-1)/2;
        % padd array via replicate
        Nsequence = padarray(Nsequence,[0 h],'replicate','both');
        Nsequence2 = [];
        for e = 1:(size(Nsequence,2)-(2*h))
            % apply QR decompose to the thread snip
            [Q,R] = qr(Nsequence(:,e:(e+2*h)));
            % if upper store - get the upper tri
            if upper
                tmp = [];
                for k = 0:size(R,1)
                    tmp = [tmp;diag(R,k)];
                end
            else
                % else get the diag
                tmp = diag(R);
            end
            % store the output sequence
            Nsequence2(:,e) = tmp(:);
        end

        threadMetric = [threadMetric;Nsequence2];
    end
    
    
end