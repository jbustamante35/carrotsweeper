function [model] = build_dataModel(data,gS,para)        
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
            % obtain subset of domain and construct model
            do_subSet = myT(data.domain(idx,:));
            do_model = do_subSet.decompose();
            do_latent = do_model.project(do_subSet);
            
            
            %%%%%
            % obtain subset of codomain and construct model
            codo_subSet = myT(data.codomain(idx,:));
            codo_model = codo_subSet.decompose();
            codo_latent = codo_model.project(codo_subSet);
            
            
            
        end
    end
end