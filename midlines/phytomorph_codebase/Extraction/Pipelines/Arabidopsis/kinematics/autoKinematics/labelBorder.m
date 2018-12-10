for e = 1:numel(SET)
    I(:,:,e) = imread(SET{e}{1});
    imshow(I(:,:,e),[]);
    drawnow
end
%% filter images
clear sI;
for e = 1:size(I,3)
    sI(:,:,e) = imfilter(double(I(:,:,e)),fspecial('gaussian',[41 41],15),'replicate');
    e
end
%% create edge image
E = zeros(size(I));
parfor e = 1:size(I,3)
    [d1 d2] = gradient(double(I(:,:,e)));
    E(:,:,e) = (d1.^2 + d2.^2).^.5;
    e
end
%% straight edge
sE = zeros(size(I));
parfor e = 1:size(I,3)
    sE(:,:,e) = edge(I(:,:,e));
    e
end
%% get border streams
PZ = 21;
stream = getBorderStream(E(:,:,1),PZ);
stream2 = zeros([size(stream) size(E,3)]);
stream = zeros([size(stream) size(E,3)]);
parfor e = 1:size(E,3)
    stream(:,:,e) = getBorderStream(E(:,:,e),PZ);
    stream2(:,:,e) = getBorderStream(sI(:,:,e),PZ);
    %stream2(:,:,e) = stream2(:,:,e) - min(min(stream2(:,:,e)));
    %stream2(:,:,e) = stream2(:,:,e) / max(max(stream2(:,:,e)));
    e
end
[streamLocs] = getBorderLocations(I(:,:,1),PZ);
%%
sz = size(stream);
stream = reshape(stream,[sz(1) sz(2)*sz(3)]);
stream2 = reshape(stream2,[sz(1) sz(2)*sz(3)]);
dC = [mean(stream);std(stream,1,1);mean(stream2);std(stream2,1,1)]';
stream = reshape(stream,sz);
stream2 = reshape(stream2,sz);
dCv = reshape(dC,[size(dC,1)/sz(3) sz(3) 4]);
dCv = permute(dCv,[1 3 2]);
%[dS dC dU dE dL dERR dLAM] = PCA_FIT_FULL_T(stream,10);
%% smooth dC
dCs = reshape(dC,[size(dC,1)/sz(3) sz(3) 4]);
dCs = permute(dCs,[1 3 2]);
for e = 1:size(dCs,3)
    dCs(:,:,e) = imfilter(dCs(:,:,e),fspecial('average',[101 1]),'circular');
end
%% make dif
for e = 1:size(dC,3)
    dCg(:,:,e) = abs(gradient(dCs(:,:,e)')');
 %   for s = 1:size(dCg,2)
 %       dCg(:,s,e) = bindVec(dCg(:,s,e))*max(dCs(:,s,e));
 %   end
end
%% make total
dCt = cat(2,dCs,dCg);
dCt = ipermute(dCt,[1 3 2]);
dCt = reshape(dCt,[size(dCt,1)*size(dCt,2) 8]);
%% transpose C
dC = dC';
%% kmeans
learnDims = [1 3 5 7];
learnDims = [1 2 3 4 5 6 7];
kidx = kmeans((dCt(:,learnDims(1))),2);
kidx = reshape(kidx,[1 sz(2) sz(3)]);

kidx = kmeans((dCt(:,learnDims)),2);
kidx = reshape(kidx,[1 sz(2) sz(3)]);

kidx2 = kmeans((dCt(:,learnDims(2))),2);
kidx2 = reshape(kidx2,[1 sz(2) sz(3)]);
close all
plot(kidx(1,:,1))
%% hkmeans
func{1}.phi = @(x,para)label_2(x,para);
func{1}.para.value = 2;
func{2}.phi = @(x,para)label_2(x,para);
func{2}.para.value = 2;
groups = ones(size(dCt,1),1);
groups = hLabel(dCt(:,learnDims),groups,func);
%% segment on patches via gmm
k = 2;
learnDims = [1 3 4];
learnDims = [1 2 3 4 5 6 7];
opt = statset('display','iter','MaxIter',500);
obj = gmdistribution.fit((dCt(:,learnDims)),k,'Options',opt);
%%
y = [];
for g = 1:size(obj.mu,1)
    y(:,g) = mvnpdf((dCt(:,learnDims)),obj.mu(g,:),obj.Sigma(:,:,g));
end
[~, kidx] = max(y,[],2);
kidx = reshape(kidx,[1 sz(2) sz(3)]);
close all
plot(kidx(1,:,10))
%% look at labels for init guess
clear X
[X(:,1) X(:,2)] = ind2sub(size(I),streamLocs);
CL = {'r.' 'g.' 'b.' 'y.'};
for e = 1:size(I,3)
    imshow(I(:,:,e),[]);
    hold on
    %plot(X(:,2),X(:,1),'go');
    for k = 1:3
        yidx = find(kidx(1,:,e) == k);
        plot(X(yidx,2),X(yidx,1),CL{k});
        hold on
    end
    waitforbuttonpress
    hold off
end
%%
kidxT = kidx;
rootGrp = 3;
backgroundGrp = 2;
borderGroup = 1;

background_node_start = hmm_node('background_start');
background_node_end = hmm_node('background_end');
border1_node = hmm_node('border1');
root_node = hmm_node('root');
border2_node = hmm_node('border2');



background_to_background = constantTransitionFunction(.90);
background_to_border1 = constantTransitionFunction(.10);

background_mean_patch = mean(dC(kidxT==backgroundGrp,learnDims));
background_cov_patch = cov(dC(kidxT==backgroundGrp,learnDims));

backgroundD = myProb(background_mean_patch,background_cov_patch);
background_node_start.attachDistribution(backgroundD,1);

background_node_start.attachNode(background_node_start,background_to_background);
background_node_start.attachNode(border1_node,background_to_border1);




stayBorder1 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
leaveBorder1 = heavisideTransitionFunction(20,@(x,y)ge(x,y));

border_mean_patch = mean(dC(kidxT==borderGroup,learnDims));
border_cov_patch = cov(dC(kidxT==borderGroup,learnDims));

borderD = myProb(border_mean_patch,border_cov_patch);
border1_node.attachDistribution(borderD,1);

border1_node.attachNode(border1_node,stayBorder1);
border1_node.attachNode(root_node,leaveBorder1);





root_to_root = constantTransitionFunction(.6);
root_to_border2 = constantTransitionFunction(.4);

root_mean_patch = mean(dC(kidxT==rootGrp,learnDims));
root_cov_patch = cov(dC(kidxT==rootGrp,learnDims));

rootD = myProb(root_mean_patch,root_cov_patch);
root_node.attachDistribution(rootD,1);

root_node.attachNode(root_node,root_to_root);
root_node.attachNode(border2_node,root_to_border2);





stayBorder2 = heavisideTransitionFunction(20,@(x,y)lt(x,y));
leaveBorder2 = heavisideTransitionFunction(20,@(x,y)ge(x,y));

borderD2 = myProb(border_mean_patch,border_cov_patch);
border2_node.attachDistribution(borderD2,1);

border2_node.attachNode(border2_node,stayBorder1);
border2_node.attachNode(background_node_end,leaveBorder1);


backgroundD2 = myProb(background_mean_patch,background_cov_patch);
background_node_end.attachDistribution(backgroundD2,1);

background_to_background_end = constantTransitionFunction(1);
background_node_end.attachNode(background_node_end,background_to_background_end);


% create hmm
hmm = my_hmm();
hmm.addNode(background_node_start);
hmm.addNode(border1_node);
hmm.addNode(root_node);
hmm.addNode(border2_node);
hmm.addNode(background_node_end);
hmm.dn = [1 1];
%%
for r = 1:10
    dCt = reshape(dC,[size(dC,1)/sz(3) sz(3) 4]);
    dCt = permute(dCt,[1 3 2]);
    clear gidx UD
    parfor e = 1:size(dCt,3)
        sig = dCt(:,learnDims,e)';
        observation_labels = [ones(size(sig,1),1)];
        gidx{e} = hmm.Viterbi(sig,observation_labels,1);
        UD{e} = [gidx{e};sig];
        e
    end

    observation_labels = [ones(3,1)];
    hmm.update(UD,observation_labels);
end
%%
close all
clear X
[X(:,1) X(:,2)] = ind2sub(size(I),streamLocs);
CL = {'r.' 'g.' 'b.' 'y.' 'r.'};
for e = 1:100%size(I,3)
    imshow(I(:,:,e),[]);
    hold on
    %plot(X(:,2),X(:,1),'go');
    for k = 1:5
        %yidx = find(kidx(1,:,e) == k);
        yidx = find(gidx{e} == k);
        plot(X(yidx,2),X(yidx,1),CL{k});
        hold on
    end
    waitforbuttonpress
    hold off
end
