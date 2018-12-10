inFilePath = '/mnt/spaldingdata/nate/communications/coldGrant/B73 vs Mo17/';  
inFilePath = '/mnt/spaldingdata/nate/communications/coldGrant/swelling/';
inFilePath = '/mnt/spaldingdata/nate/communications/coldGrant/swelling/History/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan for new images
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%% sort SET
for e = 1:numel(SET)
    N = [];
    for img = 1:numel(SET{e})
        [p n ex] = fileparts(SET{e}{img});
        N(img) = str2num(n);
    end
    [N sidx] = sort(N);
    SET{e} = SET{e}(sidx);
end
%%
I = imread(SET{1}{1});
[J mainBOX] = imcrop(I);
%%
[S A] = measureStack(SET{1},[500 500],mainBOX,0);
%%
close all
U = mean(S(1:20,:),1);
ST = std(S(1:20,:),1,1).*size(S,1).^-.5;
errorbar(((1:numel(U))*5/60)',U,ST);
ylabel(['Percent Swelling']);
csvwrite('/mnt/spaldingdata/nate/communications/coldGrant/swelling.csv',[U' ST']);
xlabel(['Time (h)']);
%% 
LEG = {'WI-Mo17-1' 'WI-Mo17-2' 'WIV-Mo17-3' 'WI-Mo17-4' 'WI-B73-1' 'WI-B73-2' 'WIV-B73-3' 'WI-B73-4'};
nS = reshape(S,[7 8 size(S,2)]);
close all
for e = 1:size(nS,2)
    tmp = squeeze(nS(:,e,:))
    U = mean(tmp,1);
    ST = std(tmp,1,1).*size(S,1).^-.5;
    errorbar(((1:numel(U))*5/60)',U,ST);
    hold all 
end
legend(LEG)
%%
close all
%G = [1 5];
%G = [2 6];
G = [3 7];
%G = [4 8];
LEG = {'Mo17' 'B73'};
nS = reshape(S,[7 8 size(S,2)]);
CL = {'b' 'r'};
close all
SNIP = 20;
DATA = [];
for e = 1:numel(G)
    tmp = squeeze(nS(:,G(e),:));
    dt = [];
    for i = 1:size(tmp,1)
        dt(i,:) = cwt(tmp(i,:),5,'gaus1');
    end
    ridx = find(any(dt(:,SNIP:end-SNIP) > .01,2));
    tmp(ridx,:) = [];
    U = mean(tmp,1);
    ST = std(tmp,1,1).*size(S,1).^-.5;
    errorbar(((1:numel(U))*5/60)',U,ST,CL{e});
    hold all
    %plot(((1:numel(U))*5/60)',tmp',CL{e});
    DATA = [DATA U' ST'];
end
legend(LEG);
csvwrite('/mnt/spaldingdata/nate/communications/coldGrant/swelling2.csv',DATA);
%% make quick movie
figure;
[J BOX] = imcrop(I,[]);
%%
close all
h = figure;
movie = axes();
graph = axes();
writerObj = VideoWriter('swelling.avi','Uncompressed AVI');
%writerObj.Quality = 100;
open(writerObj);
for e = 1:numel(SET{1})
    ROWS = round([BOX(1) (BOX(1)+BOX(3))]);
    COLS = round([BOX(2) (BOX(2)+BOX(4))]);
    tI = imread(SET{1}{e},'PixelRegion',{COLS ROWS});    
    [boundary currentArea(e)] = getKernelArea(tI,10000);
    X = linspace(0,size(tI,2),numel(SET{1}))*5/60;
    axes(movie)
    imshow(tI,[])
    hold on
    plot(boundary(:,2),boundary(:,1),'r')
    hold off
    axes(graph);
    plot(X,100*miniS,'b')
    hold on;
    plot(X(e),100*miniS(e),'ro')
    axis([0 X(end) 0 15]);
    set(graph,'Color','none');
    hold off
    xlabel('Time (h)');
    ylabel('Percent Area Increase (%)');
    
    %tI = imcrop(tI,BOX);
    writeVideo(writerObj,getframe(h));
    %imshow(tI,[])
    %M(e) = im2frame(tI);
    %drawnow
    e
    numel(SET{1})
end
close(writerObj);
%%
miniS = calcPercentSwelling(currentArea);
%%
for e = 1:50:numel(SET)
    I = imread(SET{e});    
    imshow(I,[]);
    drawnow
end
%%
plot(S')
%%
close all
G1 = [];
flag = 1;
while flag == 1;    
    I = imread(SET{1});
    [J BOX] = imcrop(I);
    
    %%
    h1 =  figure;
    h2 =  figure;
    CL = [];
    for e = 5:1:900%numel(SET)
        ROW = round([BOX(1) BOX(1)+BOX(3)]);
        COL = round([BOX(2) BOX(2)+BOX(4)]);
        I = imread(SET{1},'PixelRegion',{COL ROW});
        
        
        %I = imread(SET{e});
        %I = imcrop(I,BOX);
        
        
          
        figure(h1);
        imshow(oI,[]);
        hold on
        VEC= [];
        for i = 1:numel(B)
            plot(B{i}(:,2),B{i}(:,1));
            tmp = B{i};        
            BW = poly2mask(B{i}(:,2),B{i}(:,1),size(I,1),size(I,2));
            VEC(i) = sum(BW(:));
        end
        CL = [CL;VEC];
        hold off
        drawnow
        figure(h2)
        tmp = bsxfun(@minus,CL,CL(1,:));
        tmp = bsxfun(@times,tmp,CL(1,:).^-1);
        plot(tmp);
    end
    a = questdlg('Keep');
    if strcmp(a,'Yes');
        G1 = [G1 CL];
    end
    a = questdlg('Another');
    if strcmp(a,'Yes');
        flag = 1;
    else 
        flag = 0;
    end
end
%%
sG1 = G1;
%%
sG2 = G1;
%%
close all
G1 = sG1;
WIND = 3;
G1 = imfilter(G1,fspecial('average',[WIND 1]));
init = mean(G1(1:10,:),1);
G1 = bsxfun(@minus,G1,init);
G1 = bsxfun(@times,G1,init.^-1);
G1(1:WIND,:) = [];
G1(end-WIND:end,:) = [];
plot(G1);
%%
close all

G2 = sG2;
G2(:,3) = [];
G2 = imfilter(G2,fspecial('average',[WIND 1]));
init = mean(G2(1:10,:),1);
G2 = bsxfun(@minus,G2,init);
G2 = bsxfun(@times,G2,init.^-1);
G2(1:WIND,:) = [];
G2(end-WIND:end,:) = [];
plot(G2)
%%
close all
errorbar(3/60*(WIND:size(G1,1)+WIND-1)',mean(G1,2),std(G1,1,2)*size(G1,2).^-.5);
hold on
errorbar(3/60*(WIND:size(G1,1)+WIND-1)',mean(G2,2),std(G2,1,2)*size(G2,2).^-.5,'r');
legend({'Mo17 - N=6','B73 - N=4'});
%axis([0 270*3/60 0 .15])
title('IBM Parents - Percent Kernel Swelling');
xlabel('Time (hr)');
ylabel('Percent Swell (%)');
%%
