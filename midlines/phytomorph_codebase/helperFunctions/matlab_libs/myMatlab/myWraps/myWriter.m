function [fileList] = myWriter(kPath,data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base path
    % loop over the data structure
    % kPath is base
    % data(e).fileName is the name of the file to spool
    % data(e).d is the data to spool to the file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(data)
        % obtain the extension
        [p,n,ex] = fileparts(data(e).fileName);
        % generate file name
        fileList{e} = [kPath data(e).fileName];
        
        % NEW
        fileList{e} = data(e).fileName;
        
        % swtich on the extension
        switch ex
            case '.csv'
                % incase the file name included a path
                [p,n,ex] = fileparts(fileList{e});
                mkdir(p);
                % write file
                csvwrite(fileList{e},data(e).d);
            case '.mat'
                % incase the file name included a path
                [p,n,ex] = fileparts(fileList{e});
                mkdir(p);                
                d = data(e).d;
                save(fileList{e},'d');                
            case '.tif'
                % in case the file name included a path
                [p,n,ex] = fileparts(fileList{e});
                mkdir(p);                
                saveas(data(e).d,fileList{e});                
        end
    end
end