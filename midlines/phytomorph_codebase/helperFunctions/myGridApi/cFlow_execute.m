function [] = cFlow_execute(matFile)
    fprintf(['job called with string:' matFile '\n'])
    if isdeployed
        [p,matFile,ext] = fileparts(matFile);
    end
    matFile = ['.' filesep matFile ext];
    matFile = strrep(matFile,'"','');
    fprintf(['Loading anonymous cJob from file:' matFile '\n']);
    load(matFile,'tmpJob');
    tmpJob.localExecute();
end