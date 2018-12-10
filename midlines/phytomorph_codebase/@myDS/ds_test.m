
%T = myDS;
T = myDS('testzab.db3');
%%
mksqlite('close')
%% generate tau particle
kt = T.generateNew([],'tau');
%% put data test and set flag
k1 = T.putData(rand(10,10));
d = T.getData(k1);
w = T.setFlag(k1,'hello','world');
T.getFlag(k1,'hello');
%% put data and get data and set flag and search by flag-key
k2 = T.putData(rand(10,10));
d = T.getData(k2);
w = T.setFlag(k2,'hello','world2');
T.getFlag(k2,'hello');
nuKey = T.flagSearchByKey('hello');
nuKey = T.flagSearchByValue('world');
%%
d = T.getData(k2);
o = T.get(k1,2);
T.getFlagKeys(k1);
%% put data and get data and set flag and search by flag-key
k5 = T.putData(rand(10,10));
d = T.getData(k5);
%% put data and get data and set flag and search by flag-key
k4 = T.putData(rand(10,10));
d = T.getData(k4);
%%
k3 = T.putData(rand(10,10));
d = T.getData(k3);
%% get type
type = T.getType(kt);
type
%% create fan out
lk1 = T.upsertIndirectNonMatterLink(true,{'ef','in'},kt,k1);
%%
lk2 = T.upsertIndirectNonMatterLink(true,{'ef','in'},kt,k2);
lk3 = T.upsertIndirectMatterLink(true,{'ef','in'},kt,k3);
lk4 = T.upsertIndirectNonMatterLink(true,{'ef','in'},k4,k3);
lk5 = T.upsertIndirectMatterLink(true,{'ef','in'},k5,k3);
%% create fan in
pointerList = T.getEffluxLinkList(kt);
struct2table(pointerList)
pointerList = T.getInfluxLinkList(k3);
struct2table(pointerList)
[linkKey] = T.queryNonDirectMatterLink(kt,k3);
%%
path = T.dijkstra(kt,k1,2)
%% link tests
T.updateDirectLink(true,kt,'k');
T.updateDirectLink(false,kt,'k');
T.updateDirectLink(true,'q','k');
lk = T.upsertIndirectNonMatterLink(true,kt,k1);
r = T.isNonDirectNonMatterLinked(kt,k1);
r = T.isNonDirectNonMatterLinked(kt,k2);
lk2 = T.upsertIndirectNonMatterLink(false,kt,k1);
ndm = T.generateNew('','ndml','');
ndml = T.upsertIndirectMatterLink(true,kt,k1);
T.getType(ndml);
%% try coding up example - hard code
% shape detection-solid too
% are all animals/bears pink?
% do all bears have color?
% do all animals?
% do all things?
% infer then conclude then question
% check source target order!!!
mksqlite('close')
A = myDS();
H1 = A.putData([],'');
H2 = A.putData([],'');

Bear = A.putTypeNode('bear');
Animal = A.putTypeNode('animal');
Thing = A.putTypeNode('thing');
Pink = A.putTypeNode('pink');
Furry = A.putTypeNode('furry');
Color = A.putTypeNode('color');
Texture = A.putTypeNode('texture');

A.propLink(Bear,H1,'isa');
A.propLink(Bear,H2,'isa');
A.propLink(Animal,Bear,'typeof');
A.propLink(Animal,Thing,'typeof');
A.propLink(Pink,H1,'ism'); %
A.propLink(Color,Pink,'typeof'); 
A.propLink(Color,Bear,'hasa'); %
A.propLink(Furry,H1,'ism');
A.propLink(Texture,Furry,'typeof');
A.propLink(Texture,Bear,'hasa'); %



%{
B = myDS();
H1 = B.putData([],'');
H2 = B.putData([],'');


Bear = B.putTypeNode('bear');
Animal = B.putTypeNode('animal');
Thing = B.putTypeNode('thing');
Pink = B.putTypeNode('pink');
Furry = B.putTypeNode('furry');
Color = B.putTypeNode('color');
Texture = B.putTypeNode('texture');

B.propLink(Bear,H1,'isa');
B.propLink(Bear,H2,'isa');
B.propLink(Animal,Bear,'typeof');
B.propLink(Animal,Thing,'typeof');
B.propLink(Pink,H1,'ism'); %
B.propLink(Color,Pink,'typeof'); 
B.propLink(Color,Bear,'ism'); %all bears are pink
B.propLink(Furry,H1,'ism');
B.propLink(Texture,Furry,'typeof');
%}



A.typeExists('bear');
parentPaths = A.getParentTypes(Bear);
[type] = A.detailPath(parentPaths);
parentPaths = A.getParentTypes(H1);
[type] = A.detailPath(parentPaths);


parentPaths = A.getPotentialProperties(H1);
[type] = A.detailNodeList(parentPaths)

%%
path = T.dijkstra(Bear,{H1,H2},2);
%%
T.putReifiedTypeNode('bear','')
%%
r = T.searchTypeNetwork('Animal');
%%
mksqlite('close')
T = myDS();
imageFileNameType = T.putTypeNode('imageFileName');
fileNameType = T.putTypeNode('fileName');
referenceFrameType = T.putTypeNode('refererenceFrame');
T.propLink(fileNameType,imageFileNameType,'typeof');
T.propLink(referenceFrameType,imageFileNameType,'hasa');

parentPaths = T.getPotentialProperties(imageFileNameType);
[type] = T.detailNodeList(parentPaths)
%%

FilePath = '/mnt/tetra/nate/fixPOP/next/20180221_Rack2_Camera6/';
FileList = {};
FileExt = {'tiff'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%%

Pink = T.putTypeNode('pink');
Furry = T.putTypeNode('furry');
Color = T.putTypeNode('color');
Texture = T.putTypeNode('texture');
T.propLink(Bear,H1,'isa');
T.propLink(Bear,H2,'isa');
T.propLink(Animal,Bear,'typeof');
T.propLink(Thing,Animal,'typeof');
T.propLink(Pink,H1,'ism');
T.propLink(Color,Pink,'typeof');
T.propLink(Furry,H1,'ism');
T.propLink(Texture,Furry,'typeof');
T.typeExists('bear');




%%
fn = T.putData('/path/to/file/name.txt','fileName');
T.getType(fn);
[linkKey] = T.link(fn,k3,'hasA');
T.getType(linkKey)
%%
% forced non-direct matter 
fprintf('forced non-direct matter TEST:\n');
myDS.generateLinkingCommand(1,0,1,kt,k2);
% non-forced non-direct matter
fprintf('non-forced non-direct matter TEST:\n');
myDS.generateLinkingCommand(0,0,1,kt,k2);
% forced direct matter
fprintf('forced direct non-matter TEST:\n');
myDS.generateLinkingCommand(1,1,1,kt,k2);
% forced direct matter
fprintf('forced non-direct non-matter TEST:\n');
myDS.generateLinkingCommand(1,0,0,kt,k2);
% forced non-direct matter
fprintf('forced direct non-matter TEST:\n');
myDS.generateLinkingCommand(1,1,0,kt,k2);
%%

%%
T.iput()
%%
mksqlite('close')
%%

q = mksqlite( 'SELECT * FROM flinks' );
%%
%%
T.put(1,'data',[],[])
T.put(1,'NOTdata',[],[])
query = mksqlite( 'SELECT * FROM flinks' );
%%
T.generateNew('nu')
%%
T.generateNew('tau')
%%
T.search({'type'},{'data'})
%%
%%
q = mksqlite( 'SELECT * FROM fChain' );