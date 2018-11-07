function [] = writeData(fileName,oPath,NP,SNIP,NPK)

    oPath = strrep(oPath,'SLASH','--');
    oPath = strrep(oPath,'SPACE',' ');

    % load the data from disk
    out = loadData(fileName,SNIP,NP,NPK);
    
    % write out angle measurement
    fn = [oPath '--angle.csv'];
    csvwrite(fn,out.angle);
    % write out  length measurement
    fn = [oPath '--length.csv'];
    csvwrite(fn,out.length);
    % write out growthrate measurement
    fn = [oPath '--growthRate.csv'];
    csvwrite(fn,out.growthRate);
    
end