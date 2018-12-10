function [worked] = myiinit(username,password,varargin)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % string for writing out the environment file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        str = ['irodsHost irods.iplantcollaborative.org\n' ...
                'irodsPort 1247\n' ...
                'irodsUserName %var0%\n'...
                'irodsZone iplant\n'];
            
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
        fn = [opath '/.irodsEnv'];
        [fid msg] = fopen(fn, 'w');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % spool the environmental data file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid, strrep(str,'%var0%',iuser));
        fclose(fid);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call init irods command and pipe the password
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [statu,result] = system(['echo ''' password ''' | iinit'] );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if result contains string failed then signal fail
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = strfind(result,'failed');
        worked = isempty(fidx);
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:myiinit.m******\n']);
    end
end