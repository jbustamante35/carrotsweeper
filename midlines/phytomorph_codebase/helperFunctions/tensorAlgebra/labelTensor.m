function [P idx] = labelTensor(T,funcH_list,para_list)
    try
        %%%%%%%%%%%%%%%%%%
        % recursive call to label the domain of a tensor
        % via the set of functions in funcH_list                
        %%%%%%%%%%%%%%%%%%
        % recursive call to function list
        % init T as the input to the chain
        out{1} = T;
        for i = 1:numel(funcH_list)
            out{i+1} = funcH_list{i}(out{i},para_list{i});
        end
        % extract the co-ordinates for the selected point
        P = out{1}(out{end},:);
        idx = out{end};
    catch ME
        ME.message
    end
end