function [] = configurePeak(algoName)
    if isdeployed()
        setenv('isDeployed','true');
        setenv('phytoMorphAlgo',algoName);
        setenv('toPush','false');
    end
end