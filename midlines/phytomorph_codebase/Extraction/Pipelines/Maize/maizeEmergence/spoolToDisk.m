function [] = spoolToDisk(outFile,miniStack,miniMask)
    save(outFile,'miniStack','miniMask');
end