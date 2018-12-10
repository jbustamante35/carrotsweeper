function [] = syncMaizeDataFromIplant(user)
    % sync from iplant to local drive
    CMD = ['/mnt/scratch1/phytoM/services/maizeScannerDataSync_fromiPlant.sh ' user];
    [status,result] =  system(CMD,'-echo');
end