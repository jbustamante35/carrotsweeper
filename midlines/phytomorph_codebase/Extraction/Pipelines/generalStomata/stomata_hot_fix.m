
%% scan for maize NMS files - RILS
FilePath = '/mnt/tetra/nate/RILs/';
iFileList = {};
FileExt = {'nms'};
iFileList = gdig(FilePath,iFileList,FileExt,1);
%%
tmpI = imread(iFileList{1});
% define the radius to look over
R = [0 40];
% define the number of points
N = [(R(2)-R(1)) round(2*pi*(R(2)))];
% create the sample grid
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
% create the apply grid
[d1,d2] = ndgrid((R(2)+1):(size(tmpI,1)-R(2)),(R(2)+1):(size(tmpI,2)-R(2)));
% generate the index position(s) of the apply grid
indexPosition = sub2ind(size(tmpI),d1(:),d2(:));
% generate the sample patch
Xd = n1.*cos(n2);
Yd = n1.*sin(n2);
% store the domain
Domain = [Xd(:) Yd(:)];
% create position index list function
%plFunc = @(I)(1:numel(I));
%plFunc = @(I)indexPosition;




plFunc = @(I)generateImageDomain(I,R(2));
% generate the sampler function
samplerFunction = @(I,P)myInterp2Sampler(I,P,Domain,N);
% generate the preprocess function
preprocessFunction = @(I,P)fftPatch(I,P,samplerFunction,[1:15]);
% generate the network apply

newSZ = size(d1);
networkApply = @(I,P)applyNetworkToPoint(I,P,preprocessFunction,trainedNet_FM2,newE,newU);
%%

stomataMaster(iFileList{70},networkApply,plFunc,newSZ,'./output/','');
%%
fn = '/iplant/home/leakey_cyverse/maizeData/stomataTopoData/RILs/451 leaf2-4.nms';
func(fn)
%% publish stomataMaster
% define the radius to look over
R = [0 40];
% define the number of points
N = [(R(2)-R(1)) round(2*pi*(R(2)))];
% create the sample grid
[n1,n2] = ndgrid(linspace(R(1),R(2),N(1)),linspace(-pi,pi,N(2)));
% create the apply grid
[d1,d2] = ndgrid((R(2)+1):(size(tmpI,1)-R(2)),(R(2)+1):(size(tmpI,2)-R(2)));
% generate the index position(s) of the apply grid
indexPosition = sub2ind(size(tmpI),d1(:),d2(:));
% generate the sample patch
Xd = n1.*cos(n2);
Yd = n1.*sin(n2);
% store the domain
Domain = [Xd(:) Yd(:)];
freqToKeep = [1:15];
newSZ = size(d1);
func = @(X)stomataMaster(X,trainedNet_FM2,R,Domain,N,freqToKeep,newE,newU,newSZ,'./output/','');
pF = partialFunction(func,'maizeStomataNNapp');
pF.publish();

