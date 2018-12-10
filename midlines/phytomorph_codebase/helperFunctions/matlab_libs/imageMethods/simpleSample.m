function [S] = simpleSample(I,P,window)
    h = window(1);    
    w = window(2);
    
    % switch on class of I
    if ~isa(I,'char')
        % create sample
        S = zeros(2*h+1, 2*w+1, size(I,3), size(P,1),'uint8');
        % I is in memory
        for f = 1 : size(P,1)    
            S(:,:,:,f) = I((round(P(f,1)) - h : round(P(f,1)) + h),(round(P(f,2)) - w : round(P(f,2)) + w),:);                    
        end
    else
        info = imfinfo(I);
        S = zeros(2*h+1, 2*w+1, numel(info.BitsPerSample), size(P,1),'uint8');
        % I is in on disk
        paraPatch{2} = window(1);
        paraPatch{3} = window(2);
        for f = 1 : size(P,1)
            paraPatch{1} = P(f,:);
            S(:,:,:,f) = myReader(I,paraPatch);
        end
    end
    
end