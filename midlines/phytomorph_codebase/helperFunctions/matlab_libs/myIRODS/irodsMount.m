function [local_irodsMNT,user,remote_iRODS,remote_iRODS_return,pasword] = irodsMount(moniker)
    try 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init vars
        if nargin == 0;moniker = 'phiRODS';end
        remote_iRODS = '/iplant/home';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get home directory
        [res homePTH] = system('echo $HOME');
        homePTH(end) = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mount irods drive
        local_irodsMNT = [homePTH filesep moniker];
        [local_irodsMNT,user,pasword] = iRODSmnt(local_irodsMNT);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate irods path
        remote_iRODS = [remote_iRODS '/' user '/'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate return path for data
        remote_iRODS_return = [remote_iRODS 'phytoMorph/'];
        
    catch ME
        ME
    end
    fprintf(['done@mmount@privateIRODS\n']);
end