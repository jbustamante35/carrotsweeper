inFilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan for new images
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%%
STACK = [];
NM = [];
TARGET = [];
for e = 1:300%numel(SET)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % upper left corner for cropping out seedlings
    SEEDLING_WINDOW = [50 200 50 100];  
    I = imread(SET{e}{1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in the window which should contain the seedlings    
    I = double(I(SEEDLING_WINDOW(1):end-SEEDLING_WINDOW(3),SEEDLING_WINDOW(2):end-SEEDLING_WINDOW(4)))/255;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate cropboxes
    CROPBOX = generateCropBox(I);
    I = imresize(I,.25);
    %{
    sz = size(I);
    %tmp = im2col(I,[11 11],'sliding');
    %STACK = [STACK,tmp];
    %}
    STACK = [STACK,I(:)];
    %NM = [NM;e*ones(size(tmp,2),1)];
    NM = [NM;e];
    %answer = inputdlg('Number of objects');
    %imshow(I,[]);
    %V = impixel();
    %TARGET = [TARGET;str2num(answer{1})*ones(size(tmp,2),1)];
    %TARGET = [TARGET;size(V,1)];
    TARGET = [TARGET;numel(CROPBOX)];
    e
end

%%
STACK = double(STACK);
[XL,YL,XS,YS,BETA] = plsregress(STACK',TARGET,5);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upper left corner for cropping out seedlings
SEEDLING_WINDOW = [50 200 50 100];    
for e = 100:200
    
    
    
    I = imread(SET{e}{1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read in the window which should contain the seedlings    
    I = double(I(SEEDLING_WINDOW(1):end-SEEDLING_WINDOW(3),SEEDLING_WINDOW(2):end-SEEDLING_WINDOW(4)))/255;
    I = imresize(I,.25);
    %I = double(I);
    nk = [ones(1,1) I(:)']*BETA;
    imshow(I,[])
    title((num2str(round(nk))))
    drawnow
    pause(1)
    
end
%%
STACK = [];
NM = [];
TARGET = [];
BOX = 21;
sam = [];
CB = [];
for e = 1:10%numel(SET)
    I = imread(SET{e}{1});
    [MASK tmp] = createKernelEdgeMask(SET{e}{1});
    CB = cat(3,CB,tmp);
    [r c] = find(MASK);
    for k = 1:numel(c)
        sam = cat(3,sam,I(r(k)-BOX:r(k)+BOX,c(k)-BOX:c(k)+BOX));
    end
    e
end
szcb = size(CB);
CB = reshape(CB,[size(CB,1)*size(CB,2) size(CB,3)]);
%%
[S C U E L ERR LAM] = PCA_FIT_FULL(CB',2);
%% generate kernel mask(s)
for e = 1:10
     [MASK] = createKernelEdgeMask(SET{e}{1});
end
%% 
[S C U E L ERR LAM] = PCA_FIT_FULL(STACK',1);
%%
iterPCA(STACK',1,5)
%%
[IN] = ahe(C,TARGET',NM);
%% look at IN struct
UQ = unique(NM);
for u = 1:numel(UQ)
    tS = IN.eval(NM==UQ(u),:);
    tS = col2im(tS',[11 11],sz,'sliding');
    imshow(tS,[])
    title(sum(tS(:)))
    drawnow
    
    waitforbuttonpress
    tS = S(NM==UQ(u),:);
    tS = col2im(tS(:,(end-1)/2)',[11 11],sz,'sliding');
    imshow(tS,[])
    drawnow
    waitforbuttonpress
end
%%
func{1}.phi = @(x,para)label_2(x,para);
func{1}.para.value = 2;
func{2}.phi = @(x,para)label_2(x,para);
func{2}.para.value = 2;
func{3}.phi = @(x,para)label_2(x,para);
func{3}.para.value = 2;
func{4}.phi = @(x,para)label_2(x,para);
func{4}.para.value = 2;
groups = ones(size(STACK,2),1);
groups = hLabel(STACK',groups,func);
%% look at group labels
UQ = unique(NM);
for u = 1:numel(UQ)
    tS = groups(NM==UQ(u),:);
    tS = col2im(tS',[11 11],sz,'sliding');
    RGB = label2rgb(tS);
    imshow(RGB,[])
    drawnow
    waitforbuttonpress
end
%% junk
UQ = unique(NM);

for u = 1:numel(UQ)
    tC = C(NM==UQ(u),:);
    tS = S(NM==UQ(u),:);
    tS = col2im(tS(:,(end-1)/2)',[11 11],sz,'sliding');
    imshow(tS,[])
    drawnow
    waitforbuttonpress
    RGB = [];
    for i = 1:3%size(C,2)
        G = col2im(tC(:,i),[11 11],sz,'sliding');
        RGB = cat(3,RGB,G);
        imshow(G,[]);
        drawnow
        pause(.1)
    end
    imshow(RGB,[])
    waitforbuttonpress
end




