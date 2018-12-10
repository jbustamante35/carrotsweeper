function [results] = applyFuncToLocation(I,funcToApply,positionFunc,toDisk,verbose)
    % default to verbose
    if nargin <= 4
        verbose = true;
    end


    % read image if char
    if ischar(I)
        I = imread(I);
    end
    
    % get positions
    position = positionFunc(I);
    
    
    % pre-allocate
    if ~toDisk
        resultTmp = funcToApply(I,position(1));
        results = zeros(numel(position),numel(resultTmp));
    end
    
    % iterate
    tm = clock;
    fprintf('**********************************************\n');
    fprintf(['Starting application of function to image.\n']);
    fprintf(['Func:' func2str(funcToApply) '\n']);
    fprintf('**********************************************\n');
    parfor e = 1:numel(position)
        if verbose
            tic
            fprintf(['call start:' num2str(e) '.\n']);
        end
        tmp = funcToApply(I,position(e));
        results(e,:) = tmp(:);
        if verbose
            fprintf(['call end:' num2str(e) ':' num2str(toc) '.\n']);
        end
    end
    tm = etime(clock,tm);
    fprintf('**********************************************\n');
    fprintf(['Ending application of function to image.\n']);
    fprintf(['Func:' func2str(funcToApply) '\n']);
    fprintf(['Func in:' num2str(tm) '\n']);
    fprintf('**********************************************\n');
    
    % resize if needed
    results = reshape(results,[numel(position) size(resultTmp)]);
end

%{


%}