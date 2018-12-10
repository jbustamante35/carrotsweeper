ch1 = CH('rand',4,{1 1},{10 100});
xM1 = xM(ch1,{2 10},'Uniform','Uniform');
xM1.gen_xM();
gc = germCell();

%% main
close all
[CHp] = generateCHpara();
%% generate many xOverM for parameters
TxM = {};
MxM = {};
% generate NC number of cross over matrix
NC = 10;
for n = 1:NC
   [MxM{n}] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
end
%% single miosis event
close all
% generate a CH and the CN name
[C] = generateCH(CHp.cLv,4,'rand');
fprintf(['generated:' num2str(cN) ':ch\n']);
% miosis
[G] = miosis(MxM{1},C);
%% separate and fuse
[G0,G1] = selectGam(C);
[Cn] = fuseGam(G0,G1);
all(Cn.ch == C.ch)
%% plot check
close all
plot(C.ch,'b');
hold on
plot(C.cn(:,1),'g');
plot(.5*G.ch,'k')
plot([xVec{:}],'r')
%% self test - different ways to perform self tests
genTI = [];
close all
for g = 1:10000
    
    [CHp] = generateCHpara({1 1},{100 1000});
    close all;
    [C] = generateCH(CHp.cLv,4,'rand');
    p = [];
    [xM] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);

    for e = 1:100
        % generate two games
        [G1] = miosis(xM,C);
        [G_1{1},G_1{2}] = selectGam(G1);
        % select gam from two
        Gs1 = G_1{round(rand(1))+1};

        % generate two games
        [G2] = miosis(xM,C);
        [G_2{1},G_2{2}] = selectGam(G2);
        % select gam from two
        Gs2 = G_2{round(rand(1))+1};


        % fuse games
        [Cn] = fuseGam(Gs1,Gs2);

        % measure inbreedness
        [p(e) sites] = percentCHsim(C,Cn);
        C = Cn;

    end

    tmp = find(p == 1);
    tmp = find(p > .99);
    genTI(g,:) = [CHp.cLv(1),tmp(1)];
    plot(genTI(:,1),genTI(:,2),'r.')
    drawnow
end
plot(p)
%% generate large chrosome
[CHp] = generateCHpara({1 1},{100000 1000000});
[C] = generateCH(CHp.cLv,4,'rand');
[xM] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
% generate two games
[G1] = miosis(xM,C);
%% generate fitness vs population
[CHp] = generateCHpara({1 1},{30000 30000});
FNp = generatefitnessFunctionPara('Uniform',CHp.cLv,3);
nPOP = 100;
clear POP
for e = 1:nPOP
    POP(e) = generateCH(CHp.cLv,4,'rand');
end
POPsz = [];
close all
for t = 1:50
    POP = deathCycle(FNp,POP,.3);
    POP = birthCycle(CHp,POP,.3);
    POPsz(t) = numel(POP);
    plot(POPsz,'r')
    drawnow
end
%% test id by decs
[CHp] = generateCHpara({1 1},{10000 10000});
[C1] = generateCH(CHp.cLv,4,'rand');
[C2] = generateCH(CHp.cLv,4,'rand');
[xM] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
[Cn] = crossInv(C1,C2,CHp);
[ibdM] = idb(C1,C2);
ibdM == 0
[ibdM] = idb(Cn,C2)
%% test kinship matrix growth and deat for large rand pop
[CHp] = generateCHpara({1 1},{100 100});
FNp = generatefitnessFunctionPara('Uniform',CHp.cLv,2,1);
nPOP = 100;
clear POP
[xM] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
for e = 1:nPOP
    POP(e) = generateCH(CHp.cLv,2,'rand');
end
POPsz = numel(POP);
close all
t = 1;
while t ~= 100 & POPsz(end) < 300;
    POP = deathCycle(FNp,POP,1);
    POP = birthCycle(CHp,POP,1,xM,xM);
    t = t + 1;
    POPsz(t) = numel(POP);
    plot(POPsz,'r')
    drawnow
end
%% test kinship matrix growth and deat for large rand pop
[CHp] = generateCHpara({1 1},{100 100});
FNp = generatefitnessFunctionPara('Uniform',CHp.cLv,2,1);
nPOP = 2;
clear POP
[xM] = generateXover(CHp.xSites,CHp.cLv,CHp.exeProb);
for e = 1:nPOP
    POP(e) = generateCH(CHp.cLv,2,'inbred');
end
POPsz = numel(POP);
close all
t = 1;
while t ~= 100 & POPsz(end) < 300;
    POP = deathCycle(FNp,POP,1);
    POP = birthCycle(CHp,POP,1,xM,xM);
    t = t + 1;
    POPsz(t) = numel(POP);
    plot(POPsz,'r')
    drawnow
end
%% kinship and LD matrix
[k] = kinship(POP);
%% calc lilnkage dis
[ld] = calcLD(POP,'');
%% calc recombination frequency
sidx = randi(numel(POP),10,1);
sPOP = POP(sidx);
[md] = extractMarkerDataFromPop(sPOP);
cM = cosegregationMetric(md);
v = mean(cM(:));
for e = 1:size(cM,1)
    cM(e,e) = v;
end
imshow(-cM,[]);
%%
[rf] = recombinationFrequency(md);




%% show xover
imshow(xM,[])