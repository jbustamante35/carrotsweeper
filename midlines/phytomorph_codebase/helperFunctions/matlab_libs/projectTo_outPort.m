function [] = projectTo_outPort(outPort,data)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % configure@outPort    
    %%%%%%%%%%%%%%%%%%%%%%%%
        if 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate output filenames and data
        % whole data
        data.toDisk(1).fileName = ['whole_data.mat'];
        data.toDisk(1).d = db;
        % tip angle
        data.toDisk(2).fileName = ['tip_angle.csv'];
        data.toDisk(2).d = ta;
        % length
        data.toDisk(3).fileName = ['length.csv'];
        data.toDisk(3).d = len;
        % growth rate
        data.toDisk(4).fileName = ['growth_rate.csv'];
        data.toDisk(4).d = gr;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the index for the matlab dat
        data.toStore.IDX = 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % spool@outPort
    %%%%%%%%%%%%%%%%%%%%%%%%
    [] = spool(outPort,data)
end