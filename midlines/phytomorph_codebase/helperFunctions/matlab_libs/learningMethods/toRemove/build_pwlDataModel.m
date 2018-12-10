function [model] = build_pwlDataModel(data,gS,para)        
    %%%%%%%
    % build data model for each grouping method
    for m = 1:numel(gS)
        %%%%%
        % for each group under method m
        UQ = unique(gS(m).groups);
        for u = 1:numel(UQ)
            %%%%%
            % isolate the uth group
            idx = gS(m).groups == UQ(u);
            %%%%%
            % build model of data
            model{m}{u} = buildModel_type0(data.domain(idx,:),para);
        end
    end
end