function [ret] = pbPCA_stack_loop(fileList,objSet,tD)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % denormalize the curves
    for tm = 1:numel(objSet)
        objSet{tm}.normalize(1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over objects in set
    TOT = numel(fileList);
    for tm = 1:TOT
        fprintf(['pbPCA:' num2str(tm) ':' num2str(TOT) '\n']);
        ret{tm} = pbPCA_file_loop(fileList{tm},objSet,tD,tm);
    end
    
    for obj = 1:numel(objSet)
        g{obj} = phytoAffine();
        for tm = 1:TOT
            g{obj}.insertT(ret{tm}{obj});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize the curves
    for tm = 1:numel(objSet)
        objSet{tm}.normalize(-1);
    end
end