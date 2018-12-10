dataPath = '/iplant/home/leakey_cyverse/sorghumData/stomataTopoData/Accessions_2016';  
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
FileList = {};
FileExt = {'nms'};
 for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
 end
%%
[FileList] = issueBulkTicket(FileList);
%%
remoteOutputLocation = ['/iplant/home/leakey_cyverse/quickReturn2/'];
[remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),5*numel(FileList),'write');
%%
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    if ~isempty(strfind('EF0017_1_3',nm))
        e
        nm
        FileList{e}
    end
end
%% run local
figure;
% BAD = 240;

stomata_cnnVersion(FileList{190},convnetX,'','');
%%
func = cFlow('stomata_cnnVersion');
func.setMCRversion('v920');
for e = 1:5000
    func(FileList{e},convnet2,'./output/',remoteOutputLocation);
end

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,150,150);
%%
func = cFlow('DEwrapper');
func.setMCRversion('v840');
for e = 1:100
    func('sorghumStomata',FileList{e},true);
end

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,150,150);
%%
func = cFlow('condorDeploy_ver0');
func.setMCRversion('v840');
for e = 1:numel(FileList)
    func(FileList{e},[UW;UW2],R,[sE1,sE2],[sU1',sU2'],[Weights,Weights2],[beta1,beta2],sz,remoteOutputLocation,1);
    e
end

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
    
func.submitDag(auth,150,150);