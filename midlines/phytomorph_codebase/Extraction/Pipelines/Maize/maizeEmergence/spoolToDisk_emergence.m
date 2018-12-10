function [] = spoolToDisk_emergence(outFile,miniStack,miniMask)
    save(outFile,'miniStack','miniMask');
end