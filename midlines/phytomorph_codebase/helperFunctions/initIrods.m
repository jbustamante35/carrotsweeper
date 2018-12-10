function [] = initIrods()    
    %%%%%%%%%%%%%%%%%%%%%%%    
    % mount icommands on compute nodes
    %%%%%%%%%%%%%%%%%%%%%%%
    if isdeployed
        fprintf(['Initializing to irods@iplant \n']);
        [worked] = myiinit('phytomorphuser','phyt0M0rph',1);
        if worked
            fprintf(['Initialization to irods@iPlant worked.\n']);
        else
            fprintf(['Initialization to irods@iPlant did not work.\n']);
        end
    end
end