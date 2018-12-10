%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lite test my tensor class
sz = [5,5,6];
T = myT(rand(sz));
fprintf(['************************\n']);
nsz = circshift([1 2 3],-1);
% next test - permute and invert permute
T.permute(nsz);
osz = size(T);
T.ipermute(nsz);
sz = size(T);
fprintf(['report from permute and ipermute dims:' num2str(all(osz == sz)) '\n']);
fprintf(['************************\n']);
fprintf(['************************\n']);
% next test - permute and invert
nsz = circshift([1 2 3],-1);
T.permute(nsz);
osz = size(T);
T.i();
sz = size(T);
fprintf(['report from permute and invert:' num2str(all(osz == sz)) '\n']);
fprintf(['************************\n']);
% to bson
oT = T.toBson();
nT = myT.fromBson(oT);
mT = myT(oT);
fprintf(['Tested permute\n']);
% test plus
v = myT(rand(1,4));
w = myT(rand(1,4));
z = v + w;
%% test myBasis class
BV = phytoAbasis(ones(4,4));
BV(1) = rand(4,1);
BV(2) = rand(4,1);
%% test fibreBundle
fprintf(['*****************************************************\n']);
fprintf(['Testing notation for fibreBundle:START\n']);
% make fibre and base
d = rand(2,3,4,5);
% create rank-2 base and rank-2 tensor fibre
f = myTb(2,2,d);
% index the base
f1 = f(2,3);
% assign the fibre
f(1,3) = 2;
% test the colon operation for the tensor
f(:,1);
% test from the first to the end
f(1:end,1);
% test the 1:N
sz = [4 5 2 4];
SZ = prod(sz);
H = myTb(1,4,reshape(1:SZ,sz));
% toBson
of = f.toBson();
kf = myTb(of);
fprintf(['Testing notation for fibreBundle:END\n\n']);
%% test fibreBundle(s) - rank--n - not sure about using thie structure
fprintf(['*****************************************************\n']);
fprintf(['Testing notation for fibreBundle n-Rank:START\n']);
for loop = 1:10
    for baseR = 0:3        
        for fibreR = 0:5
            fprintf(['Testing notation for baseRank:fibreRank' num2str(baseR) ':' num2str(fibreR) '...\n']);
            fprintf([num2str(baseR) ':' num2str(fibreR) '...']);
            % create sentence
            A = {'1','2','3','4','5','6','7','8','9'};
            D = ',';
            L = baseR + fibreR;
            s = genString(A,L,D);
            
            % create : sequece
            A = {':'};
            D = ',';            
            sF = genString(A,fibreR,D);
            
            % base indexd
            indS = s(1:(baseR + (baseR-1)));
            fibreSZ = s(baseR + (baseR-1)+2:end);
            fSZ = str2num(fibreSZ);
            
            % create base and fibre
            cmd =[ 'd = rand(' s ');'];
            eval(cmd);

            % create base and fibre(s) object
            f = myTbS(baseR,fibreR);
            % set the base fibre :: aka the base
            f{1}.setData(d);
            
            
            
            cmd = ['g = f{1}(' indS ');'];
            eval(cmd);                        
        end
        fprintf(['\n']);
    end
end
fprintf(['Testing notation for fibreBundle n-Rank:END\n\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% affine point test
fprintf(['*****************************************************\n']);
fprintf(['Testing notation for affinePoint:START\n']);
Ps = phytoApoint([rand(2,1);1]');
oPs = Ps.toBson();
kPs = phytoApoint(oPs);
delta = phytoApoint([rand(2,1);1]');
% view the point
h = figure;hold on;
vProps.Marker = 'o';
vProps.Color = 'g';
Ps.view(h,vProps);
Ps = Ps + delta;
vProps.Marker = 'o';
vProps.Color = 'r';
Ps.view(h,vProps);
Ps = Ps - delta;
fprintf(['Testing notation for affinePoint:END\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% affine curve test
% note that curves are stored as along fibre is along dim one
fprintf(['*****************************************************\n']);
fprintf(['Testing notation for affineCurve:START\n']);
T = linspace(0,3,100);
func{1} = @(T) [T.^2;T;ones(size(T))];
func{2} = @(T) [T.^1 + T.^2;T;ones(size(T))];
X = func{1}(T);
delta = func{2}(T);
C = phytoAcurve(X');
oC = C.toBson();
kC = phytoAcurve(oC);
D = phytoAcurve(delta');
h = figure;hold on;
C.view(h,[]);
C = C + D;
C.view(h,[]);
fprintf(['Testing notation for affineCurve:END\n']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% affine transformation test
D = eye(3);
D(1:2,end) = [3 4];
C = phytoAaffine(D);
h = figure;hold on;
C.view(h,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% affine domain
para.type = 'disk';
para.value{1} = [0 20 200];
para.value{2} = [-pi pi 200];
D = phytoAdomain(para);
D.generateDomain();
h = figure;hold on;
D.view(h,[],[]);
% test displacement via mult
delta = eye(3);
delta(1:2,end) = [3 4];
delta = phytoAaffine(delta);
D = delta*D;
D.view(h,[],[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test with cloud of affine points
NP = 30;
Ps = myHS_X('phytoApoint');
for e = 1:NP
    p = phytoApoint([rand(2,1);1]');
    Ps{e} = (p);
end
h = figure;hold on;
Ps.view(h,[],[]);

vProps.Marker = 'o';
vProps.Color = 'g';
Ps.view(h,[],vProps);

vProps.Marker = '*';
vProps.Color = 'b';
Ps.view(h,[],vProps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% graph
G = myG();
% generate X rand nodes
N = 100;
for e = 1:N
    tN = myN();
    pt = phytoApoint([rand(2,1);1]');        
    tN{1} = pt;
    G.putNode(tN);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;hold on;
vProps.nodes = [];
vProps.edges = [];
G.view(h,[],vProps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simple test of graph with nodes having only points
nodeType = 'points';
%nodeType = 'affine';
mag = 10;
%%%%%%%%%%%%%%%%
% graph
G = myG();

%%%%%%%%%%%%%%%%
% generate X rand node
N = 100;
for e = 1:N
    
    tN = myN();
    
    switch nodeType
        case 'points'
            
            v1 = [rand(2,1);1]';
            pt = phytoApoint(v1);
            
        case 'affine'
            E = eye(3);
            
            v1 = rand(2,1)-.5;
            v2 = rand(2,2)-.5;
            
            E(1:2,3) = mag*v1;
            E(1:2,1:2) = mag*v2;
            pt = phytoAaffine(E);
            
    end
    tN{1} = pt;
    G.putNode(tN);
end


%%%%%%%%%%%%%%%%
% generate M connections
M = 10;
for e = 1:M
    r1 = round(1 + (N-1).*rand(1));
    r2 = round(1 + (N-1).*rand(1));
    tE = myEdge();
    tE.setSource(G.getNode(r1));
    tE.setTarget(G.getNode(r2));
    G.putEdge(tE);
    SL(e) = r1;
end

h = figure;hold on;
G.view(h,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view above in motion
h = figure;

p = @(t).2*[cos(t),sin(t),0];
TH = linspace(-pi,pi,30);

for e = 1:numel(TH)
    
    dV = phytoApoint(p(TH(e)));
    
    % get source node 1 SL contains the list of source nodes
    n = G.getNode(SL(1));
    
    
    pt = n{1};
    pt = pt + dV;
    n{1} = pt;
    G.putNode(n,SL(1));
    hold on
    G.view(h,[]);
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test sampler code for point
p_disk.value{1} = [0 200 200];
p_disk.value{2} = [-pi pi 200];
p_disk.type = 'disk';

p_box.value{1} = [-100 100 200];
p_box.value{2} = [-400 400 400];
p_box.type = 'box';

disk = phytoAdomain(p_disk);
box = phytoAdomain(p_box);

disk.generateDomain();
box.generateDomain();

%%%%%%%%% read image and gather point
iPth = '/mnt/scratch1/deploy/sampleData/shortStack/phytoMorph/rootAnalysis/graviTropism_kinematics/071811 adgspl2-10/000525.TIF';

I = myReader(iPth);
[r c v] = impixel(I);
close all;
for e = 1:numel(r)
    pt = phytoApoint([r;c;ones(size(1))]');
    P = pt.sample(iPth,box);
    P.view(h,[],[]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test sampler code for edge
p_disk.value{1} = [0 200 200];
p_disk.value{2} = [-pi pi 200];
p_disk.type = 'disk';

p_box.value{1} = [300 -300 200];
p_box.value{2} = [1 -1 400];
p_box.type = 'box';

disk = phytoDomain(p_disk);
box = phytoDomain(p_box);

disk.generateDomain();
box.generateDomain();
%%%%%%%%% read image and gather point
iPth = '/mnt/scratch1/deploy/sampleData/shortStack/phytoMorph/rootAnalysis/graviTropism_kinematics/071811 adgspl2-10/000525.TIF';
%%%%%%%%%%%%%%%%%%
I = myReader(iPth);
data = myDraw(I);
data = data{1};
% construct graph
G = myG();
for e = 1:(size(data,1)-1)
    %%% construct source
    source = myN();
    pt = phytoApoint([data(e,:) 1]);    
    source{1} = pt;
    %%% construct target
    target = myN();
    pt = phytoApoint([data(e+1,:) 1]);   
    target{1} = pt;
    %%% construct edge
    eg = myEdge();
    eg.setSource(source);
    eg.setTarget(target);
    %%% insert data into graph
    G.putNode(source);
    G.putNode(target);
    G.putEdge(eg);
end
%%%%%%%%%%%%%%%%%%
close all;
h = imshow(I,[]);
hold on;
%plot(data(:,1),data(:,2),'r*');
%%%%%%%%%%%%%%%%%%
G.view(h,[]);
graphDomain.edge = box;
graphDomain.node = disk;
G.sampleViewer(graphDomain,h);
%%%%%%%%%%%%%%%%%%
sG = G.sample(iPth,graphDomain);
sG.view(figure,[],[])
%%
tC = myHS_X('char');
tC{1} = 'a';