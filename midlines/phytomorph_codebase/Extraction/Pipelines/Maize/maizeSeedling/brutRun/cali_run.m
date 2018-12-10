dataPath = ['/iplant/home/canibas/maizeData/seedlingData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
%%
caliList = {};
FileExt = {'nef'}
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        caliList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%%
cidx = [];
for e = 1:numel(caliList)
    sidx = (strfind(caliList{e},'Jan'));
    if ~isempty(sidx)
        dayNum = str2num(caliList{e}((sidx(1)-3):(sidx(1)-2)))
        if dayNum >= 7 & dayNum <= 30
            cidx = [cidx e];
        end
    end
end
%%
tic

parfor e = 1:5
    
    fileName = caliList{cidx(numel(e))};
    smartMain_v4(fileName,'./output/',[],GMModels{3},level2,levelA,levelD,kE,kU);
   
end
 toc
 %%
 
fileName = caliList{cidx(numel(e))};
caliList = caliList(cidx);
[caliList] = issueBulkTicket(caliList);
%%
remoteOutputLocation = '/iplant/home/canibas/fastOSG/';
CMD = ['imkdir -p ' remoteOutputLocation];
system(CMD);
    
    numJobs = numel(caliList);
    [remoteOutputLocation iticket] = issueTicket(remoteOutputLocation(1:end-1),10*numJobs,'write');
   %%
caliFunc = cFlow('smartMain_v4');
caliFunc.setMCRversion('v930');
for e = 1:numel(caliList)
    caliFunc(caliList{e},'./output/',remoteOutputLocation,GMModels{3},level2,levelA,levelD,kE,kU);
end
  auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    caliFunc.submitDag(auth,50,50);
    %%
    massDownload('/iplant/home/canibas/fastOSG/', '.json','/home/nate/Downloads/caliOSG/');
    %%
    