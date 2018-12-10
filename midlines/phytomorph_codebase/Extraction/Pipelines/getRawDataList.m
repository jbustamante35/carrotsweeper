function [r] = getRawDataList(user,plantType,tissueType)
    dataPath = ['/iplant/home/' user '/' plantType 'Data/' tissueType 'Data%'];
    CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
    [o,r] = system(CMD);
    [r] = parseRecords(r);
end




%{
    getRawDataList('garf0012','maize','ear')
%}