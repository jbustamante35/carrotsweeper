function [IN] = ahe(data,target,groups)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data
    sflag = 0;
    approach = 'lite';    
    %approach = 'loop'; 
    % set the optimization parameters
    MAXcnt = 5;
    POP = 100;
    GEN = 500;
    thresh = .008;
    switch approach
        case 'lite'
            cnt = 1;
            while ~sflag
                % generate the bounds of the data
                bounds = modelBounds(data);
                options = psooptimset('PopInitRange',bounds,'PopulationSize',POP,'Display','iter','Generations',GEN,'Vectorized','on','CognitiveAttraction',2,'SocialAttraction',.1);    
                % run the optimization
                [mo{cnt},fval{cnt},exitflag,output,population,scores] = pso(@(m)ah3(m,data,target,groups),size(bounds,2),[],[],[],[],[],[],[],options);
                sflag = fval{cnt} < thresh;
                if cnt >= MAXcnt;sflag=1;end
                cnt = cnt + 1;
            end
        case 'loop'
            for e = 1:MAXcnt
                % generate the bounds of the data
                bounds = modelBounds(data);
                options = psooptimset('PopInitRange',bounds,'PopulationSize',POP,'Display','iter','Generations',GEN,'Vectorized','on','CognitiveAttraction',2,'SocialAttraction',.1);    
                % run the optimization
                [mo{e},fval{e},exitflag,output,population,scores] = pso(@(m)ah3(m,data,target,groups),size(bounds,2),[],[],[],[],[],[],[],options);
            end
        case 'slam'
                POP = 1000;
                GEN = 500;
                bounds = modelBounds(data');
                options = psooptimset('PopInitRange',bounds,'PopulationSize',POP,'Display','iter','Generations',GEN,'Vectorized','on','CognitiveAttraction',1.5,'SocialAttraction',1.5);    

                % run the optimization
                [mo,fval,exitflag,output,population,scores] = pso(@(m)ah3(m,Cdata,target,groups),size(bounds,2),[],[],[],[],[],[],[],options);
    end
    
    
    mo = cell2mat(mo');
    fval = cell2mat(fval);
    [fval sidx] = min(fval);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eval data
    [p] = gprobData(data,mo(sidx,:));
    IN.D = [];
    IN.model = mo(sidx,:);
    IN.eval = p;
    IN.fval = fval;
end

