 %%% for sampling curves
T = 30;
R = 50;
nhSZ = [T R];
[THETA RAD] = ndgrid(linspace(-pi,pi,T),linspace(0,30,R));
NH = [RAD(:).*cos(THETA(:)) RAD(:).*sin(THETA(:))]';
segmentSize = 61;    
%% load atomicFragments data for tip isoloation

FilePath = '/mnt/scratch5/maizeContours/';
FileList = {};
FileExt = {'mat'};
verbose = 0;
[FileList] = gdig(FilePath,FileList,FileExt,verbose);


D = {[]};
G = [];
segmentSize = 61;
parfor i = 1:200
    a = load(FileList{i});
    atomicFragments = a.obj.atomizeSingleContour(segmentSize,1);
    tmp = [];
    for e = 1:numel(atomicFragments)
        tmp = [tmp;atomicFragments(e).Osegment(:)'];
    end
    D{i} = tmp';
    %D{i+1} = [];
    %T(i) = 1;
    %G = [G;i*ones(numel(atomicFragments),1)];
    i
end
D = cell2mat(D)';
%% analysis of tip loading data
[S C U E L ERR LAM] = PCA_FIT_FULL(D,5);
%% get midlines
sPath = '/mnt/scratch5/maizeContours_withmidlines/';
mkdir(sPath)
parfor i = 1:numel(FileList)
    try
        [p n ex] = fileparts(FileList{i});
        nName = [sPath n ex];
        if ~exist(nName)
            root = load(FileList{i},'obj');
            root = root.obj;
            root.flipTraceDirectionToClockWise();
            root.traceMidline(segmentSize,E,U);        
            root.persist(nName);
        end
    catch ME
         fprintf(['Error in:' num2str(i) '\n']);
    end
    fprintf(['Done with:' num2str(i) '\n']);
end
%%
try
names = J(e).getKernelFileNames(sPath);
E = [];
for n = 1:numel(names)
E(n) = exist(names{n});
end
if ~all(E)
cS = J(e).extractContourSequences();
for c = 1:numel(cS)
    cS{c}.imageStack.flush();
    cS{c}.persist(names{c});
end
end

catch ME
ME
fprintf(['error on' num2str(e) '\n']);
end
%%
L = [];A = [];
for e = 1:numel(root)
    L(e,:) = root{e}.obj.getMidlineLength();
    A(e,:) = root{e}.getTipAngle();
end
%%
tmp = [];
%40
for i = 93
    a = load(FileList{i});
    %a.obj.traceMidline(segmentSize,E,U);
    a
    atomicFragments = a.obj.atomizeSingleContour(segmentSize,1);
    tmp = [];
    for e = 1:numel(atomicFragments)
        tmp = [tmp;atomicFragments(e).Osegment(:)'];
    end
    Jtmp = a.obj.contourSequence{1}.getIntersectionProfile();
    Jtmp = sum(Jtmp.*Jtmp,1).^.5;
    d1 = cwt(Jtmp,[30],'gaus1');
    d1 = abs(d1);
    Jtmp = closedCurve.windowData(Jtmp,31);
    Jtmp = reshape(Jtmp,[size(Jtmp,1)*size(Jtmp,2) size(Jtmp,3)]);
    [tC] = PCA_REPROJ(tmp,E,U);
    %[jC] = PCA_REPROJ(Jtmp',jE,jU);
    
    

    %% analysis of tip in the first
    sig = tC(:,1).*(a.obj.contourSequence{1}.centerVec(1,:)  < 0)';
    [J fidx] = max(sig);
    %{
    %% analysis of junction in the 4th
    COMP = 2;
    jC(:,COMP) = abs(jC(:,COMP));
    jidx = nonmaxsuppts(jC(:,COMP), 10);
    jidx = jidx .* (jC(:,COMP) > 0) .* (jC(:,COMP) > 10);
    jidx = find(jidx);
    % remove the tip from contending
    [~,ridx] = min(abs(jidx - fidx));
    jidx(ridx) = [];
    plot(jC(:,COMP));
    hold on;
    plot(jidx,jC(jidx,COMP),'r*');
    pidx = find((jidx - fidx) > 0);
    nidx = find((jidx - fidx) < 0);
    pdist = abs(jidx(pidx) - fidx);
    ndist = abs(jidx(nidx) - fidx);
    [~,spidx] = sort(pdist);
    [~,snidx] = sort(ndist);
    jidx = jidx([nidx(snidx(1)) pidx(spidx(1))]);
    %}
    
    %% analysis of der
    jidx = nonmaxsuppts(d1, 10);
    jidx = jidx .* (d1 > 100);
    jidx = find(jidx);
    plot(d1);
    hold on;
    plot(jidx,d1(jidx),'r*');
    pidx = find((jidx - fidx) > 0);
    nidx = find((jidx - fidx) < 0);
    pdist = abs(jidx(pidx) - fidx);
    ndist = abs(jidx(nidx) - fidx);
    [~,spidx] = sort(pdist);
    [~,snidx] = sort(ndist);
    jidx = jidx([nidx(snidx(1)) pidx(spidx(1))]);
    
    %%
    a.obj.contourSequence{1}.tagTip(fidx);
    a.obj.contourSequence{1}.tagJunctions(jidx);
    %
    close all
    a.obj.contourSequence{1}.plot(1);
    axis equal;
    hold on;
 
    
    
end
[J sidx] = min(jC(:,2));
%% inspect junction components
close all
plot(jC);
%plot(jC(:,4),'k');
%% OLD
close all
tProb = gprobData(tC,IN.model);
[J fidx] = max(tProb);
%% analysis of tip in the first
sig = tC(:,1).*(a.obj.contourSequence{1}.centerVec(1,:)  < 0)';
[J fidx] = max(sig);
%% analysis of junction in the 4th
COMP = 2;

jidx = nonmaxsuppts(jC(:,COMP), 30);
jidx = jidx .* (jC(:,COMP) > 0) .* (jC(:,COMP) > 20);
jidx = find(jidx);
% remove the tip from contending
[~,ridx] = min(abs(jidx - fidx));
jidx(ridx) = [];
plot(jC(:,COMP));
hold on;
plot(jidx,jC(jidx,COMP),'r*');
pidx = find((jidx - fidx) > 0);
nidx = find((jidx - fidx) < 0);
pdist = abs(jidx(pidx) - fidx);
ndist = abs(jidx(nidx) - fidx);
[~,spidx] = sort(pdist);
[~,snidx] = sort(ndist);
jidx = jidx([nidx(snidx(1)) pidx(spidx(1))]);
%% 
a.obj.contourSequence{1}.tagTip(fidx);
%%
close all
a.obj.contourSequence{1}.plot(1);
axis equal;
hold on;
plot(a.obj.contourSequence{1}.segment(1,fidx),a.obj.contourSequence{1}.segment(2,fidx),'k*');
plot(a.obj.contourSequence{1}.segment(1,jidx),a.obj.contourSequence{1}.segment(2,jidx),'g*');
%%
plot(tC);