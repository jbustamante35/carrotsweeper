function [worked] = myiinit(username,password,varargin)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % string for writing out the environment file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %str = ['irodsHost trixie.iplantcollaborative.org\n' ...
        str = ['irodsHost irods.iplantcollaborative.org\n' ...
                'irodsPort 1247\n' ...
                'irodsUserName %var0%\n'...
                'irodsZone iplant\n'];
        str = ['{ \n '...
                '\t"irods_host": "data.iplantcollaborative.org",\n'...
                '"irods_zone_name": "iplant",\n'...
                '"irods_port": 1247,\n'...
                '"irods_user_name": "%var0%",\n'...
                '"irods_authentication_file": "%var1%"\n'...
                '}'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out the file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init the user name
        iuser = username;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the .irods folder for storing the password and env files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin == 2
            [res opath] = system('echo $HOME');
            opath(end) = [];
            opath = [opath '/.irods'];
            mkdir(opath);
        end
        if nargin > 2
            [res opath] = system('echo $PWD'); 
            opath(end) = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make and open the environmental data file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fn = [opath '/.irodsEnv'];
        fprintf(['irods_icommands_dir:' opath '\n']);
        fn = [opath '/irods_environment.json'];
        [fid msg] = fopen(fn, 'w');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % spool the environmental data file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        newENVString = strrep(strrep(str,'%var0%',iuser),'%var1%',[opath filesep 'pwfile']);        
        fprintf(fid,newENVString);
        fclose(fid);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set var to file name
        CMD = ['export IRODS_ENVIRONMENT_FILE=' opath '/irods_environment.json'];        
        system(CMD,'-echo');
        CMD = ['export ENVIRONMENT_VAR_HOME=' opath];
        system(CMD,'-echo');        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call init irods command and pipe the password
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[status,result] = system(['echo ''' password ''' | iinit'],'-echo');
        CMD = ['IRODS_ENVIRONMENT_FILE=' opath '/irods_environment.json iinit ' password];        
        [status,result] = system(CMD,'-echo');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if result contains string failed then signal fail
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = strfind(result,'failed');
        worked = isempty(fidx);
        %CMD = ['ils -h'];
        %system(CMD);
    catch ME
        getReport(ME)
        
    end
end