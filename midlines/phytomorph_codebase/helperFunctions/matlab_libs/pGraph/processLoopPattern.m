function [] = processLoopPattern(S,loopPattern)
    idList = S.joint(loopPattern(1).anchorNode,loopPattern(1).nodeType);
    for e = 1:numel(idList)
        tmpL = loopPattern(2:end);
        tmpL(1).anchorNode = idList(e);
        processLoopPattern(tmpL);
    end
end

%{
rng('shuffle');
baseP = 'X:\nate\junk\';
L = 64;
alphabet = [48:57 65:90 97:122];
tmpP = char(alphabet(floor(length(alphabet)*rand(1,L)) + 1));
dbP = [baseP tmpP filesep];
nm = 'db';
S = store(nm,dbP);

head = S.createNode();

for i = 1:10
    nn = S.createNode();
    nl =S.makeLink(nn,head);
    nl.setProp('level','1');
end
%}