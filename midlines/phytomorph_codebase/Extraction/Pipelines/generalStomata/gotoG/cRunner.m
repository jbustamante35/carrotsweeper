function [] = cRunner(A,B,C,D)
    webPath = './';
    options = weboptions('Timeout',60);
    fprintf(['start:pulling applyAllLayers.mat\n']);
    websave([webPath 'sorghumNNapp.mat'],'https://de.cyverse.org/dl/d/5E60D0E9-2F48-4540-8065-A17EDFCF0EA3/sorghumStomataApply.mat',options);
    load([webPath 'sorghumNNapp.mat']);
    funcToCall = obj.func;
    delete([webPath 'applyAllLayers.mat']);
    funcToCall(A,B,C,D);
end