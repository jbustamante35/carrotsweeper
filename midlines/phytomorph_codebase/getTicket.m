function [ticket] = getTicket(user,type,kFile,oPath)
    tmpFile = [oPath 'lwe'];
    websave(tmpFile,'https://de.cyverse.org/dl/d/BB7687DF-6A4A-4CE4-9556-7360B8E33EB2/e');
    oFile = [oPath 'ue'];
    cmd = ['openssl rsautl -decrypt -inkey ' kFile ' -in ' tmpFile ' -out ' oFile];
    [status,cmdout] = system(cmd);
    kv = loadjson(oFile);
    ticket = kv.([user type]);
    delete(oFile);
    delete(tmpFile);
end