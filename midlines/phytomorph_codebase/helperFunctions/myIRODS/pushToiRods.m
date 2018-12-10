function [] = pushToiRods(location,fileList)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % push to irods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(location)
        % clip off ticket
        [pth,ticket,ext] = fileparts(location);
        fidx = strfind(location,'#');
        pth = location(1:(fidx(end-1)-1));
        ticket = location((fidx(end-1)+1):(fidx(end)-1));
        fprintf(['**************************************************************************\n'])
        fprintf(['Start push to irods:' num2str(numel(fileList)) ' files.\n']);tic
        fprintf(['**************************************************************************\n'])
        for e = 1:numel(fileList)
            fprintf(['*************************************\n'])
            fprintf(['Start push to irods:' num2str(e) ':' num2str(numel(fileList)) ' files.\n']);tic
            fprintf(['*************************************\n'])
            [~,tn,te] = fileparts(fileList{e});
            targetFile = [pth filesep tn te];
            targetHTTP = [targetFile '#' ticket '#'];
            targetHTTP = xform2URL({targetHTTP});
            targetHTTP = targetHTTP{1};
            cmd = ['iput -f -V -t ' ticket ' "' fileList{e} '" "' targetFile '"'];
            %cmd = ['curl -X POST ' targetHTTP ' -F uploadFile=@' fileList{e}];
            fprintf(['Push Command is:' cmd '\n']);
            [r,o] = system(cmd,'-echo');
            fprintf(['\npushing to file to irods:' tn '\n']);
            fprintf(['*************************************\n'])
            fprintf(['End push to irods:' num2str(e) ':' num2str(numel(fileList)) ' files.\n']);tic
            fprintf(['*************************************\n'])
        end
        fprintf(['**************************************************************************\n'])
        fprintf(['End push to irods:' num2str(toc) ' seconds.\n']);
        fprintf(['**************************************************************************\n'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % push to irods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end