function [data_model] = eval_pwlDataModel(data,model,para)
    % for each grouping method
    for meth = 1:numel(model)
        for u = 1:numel(model{meth})
            %%% NOTE router could be here
            % build model of data
            data_model{meth}{u} = evalModel_type0(model{meth}{u},data.domain);            
        end
    end
end