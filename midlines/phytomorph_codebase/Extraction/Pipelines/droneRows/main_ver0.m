I = imread('/home/nate/Downloads/p1_mica_20160823_index_ndvi.tif');
I = imread('/mnt/spaldingdata/Drone_Imagery/Arlington_2018/2018_07_23/2018_07_23.tif');
%%
I = I(:,:,1:3);
%%
G = rgb2gray(I);
%%
imshow(G,[]);
%%
subRange = imcrop(G);
close all
imshow(subRange,[]);
%% 
close all
sig = -mean(subRange,2);
figure;
plot(sig,'b')
sig = imfilter(sig,ones(100,1)/100,'symmetric');
hold on;
plot(sig,'r')
rangeD = imdilate(sig,ones(800,1)) == sig & sig > -80;
plot(-100*rangeD,'k')
%%
%axis([0 size(rangeD,1) 0 2]);
%{
IM = imreconstruct(sig-100,sig);
sig = sig - IM;
rangeD = sig > 0;
%}
raIDX = find(rangeD);
rangeD = repmat(rangeD>0,[1 size(subRange,2)]);
rangeD = imdilate(rangeD,strel('disk',11,0));
out = flattenMaskOverlay(subRange,rangeD);
figure;
imshow(rangeD,[]);
figure;
imshow(out,[]);
imwrite(out,'/mnt/tetra/nate/for_IOWA_CORN/large.tiff');
%%
subI = double(subI)/255;
%%
%%
close all
sig = fft(subI-mean(subI,2),[],2);
sig2 = mean(abs(sig),1);
ssig2 = imfilter(sig2,fspecial('average',[1 21]),'replicate');
figure;
plot(ssig2,'k')
hold on
plot(sig2);
figure;
sig3 = sig2 - ssig2;
plot(sig3);
sig4 = imdilate(sig3,strel('disk',51,0)) == sig3;
cutoff = 50;
sig4(1:cutoff) = 0;
fidx = find(sig4);
hold on
plot(sig4*30);
%%
BW = roipoly(I);
%%

close all
subI = imcrop(I);
subI = bindVec(subI);
imshow(double(subI),[]);
%%
%%
sig = [];
for r = 1:(numel(raIDX)-1)
    subI = double(subRange(raIDX(r):raIDX(r+1),:))/255;
    close all
    [f,phase,IDXf,sig(r,:)] = findMajorFrequency(subI,50,140^-1,[],0,[],[]);
    
end
F = mean(sig,1);
FF = std(sig,1,1);
%%

%%
for r = 1:numel(raIDX-1)
    subI = double(subRange(raIDX(r):raIDX(r+1),:))/255;
    close all
    [f,phase,IDXf] = findMajorFrequency(subI,50,(.5*140)^-1,[],0,F,FF);
    %[f,phase] = findMajorFrequency(subI,50,50^-1,0);
    [mask] = generateRowMask(subI,f,phase);
    out = flattenMaskOverlay(double(bindVec(subI)),mask>.8);
    mask2 = repmat(mean(mask,1),[size(mask,1) 1]);
    out2 = flattenMaskOverlay(double(bindVec(subI)),bindVec(mask2)<.02);
    imshow(out,[])
    figure;
    imshow(out2,[]);
    imwrite(out2,['/mnt/tetra/nate/for_IOWA_CORN/small' num2str(r) '.tiff']);
    title(num2str(f^-1))
    figure;
    plot(mean(subI,1));
    %waitforbuttonpress
end
%%
close all
blockNUM = 10;
WIN = round((1/f)*blockNUM);
BIG = zeros(size(subI));
for s = 1:(size(subI,2)-WIN)
    tempy = subI(:,s:(s+WIN));
    IDXf = round(size(tempy,2)*f)+1;
    
    [Tf,Tphase(s,:)] = findMajorFrequency(tempy,50,60^-1,IDXf,0);

    [Tmask] = generateRowMask(tempy,f,Tphase(s,:));
    Tmask = repmat(mean(Tmask,1),[size(Tmask,1) 1]);
    
    BLOCK = BIG(:,s:(s+WIN));
    BLOCK = mean(cat(3,BLOCK,Tmask),3);
    BIG(:,s:(s+WIN)) = BLOCK;
    imshow(BIG,[]);
    drawnow
    
    
    %{
    imshow(tempy,[]);
    drawnow
    %}
end
%%
close all
for e = 1:1
    

    out2 = flattenMaskOverlay(double(bindVec(subI)),bindVec(mask2)<.02);
    imshow(out2,[]);
    hold on
    sig = -mean(BIG,1);
    IM = imreconstruct(sig-.1,sig);
    sig = sig - IM;
    sig = repmat(sig > 0,[size(BIG,1) 1]);
    waitforbuttonpress
    out2 = flattenMaskOverlay(out2,sig,.5,'b');
    imshow(out2,[]);
    drawnow
    %plot(mean(BIG,1)*500+500,'r')
    waitforbuttonpress
    
    %plot(mean(mask2,1)*500+500,'b')
end
%%
ridx = find(any(BW,2));
PAD = 51;
phun = zeros(size(BW));
for r = PAD:(numel(ridx)-PAD)
    
    
    tic
    strip = I((ridx(r)-PAD):(ridx(r)+PAD),:);
    Mstrip = BW((ridx(r)-PAD):(ridx(r)-PAD),:);
    midx = find(all(Mstrip,1));
    subI = strip(:,midx);
    toc
    
    
    
    [~,BOX] = imcrop(I);
  
    
    tic
    [block] = getAblock(I,BOX);
    toc
    
    block = imfilter(block,fspecial('gaussian',[15 15],2),'replicate');
    
    tic
    [f,phase] = findMajorFrequency(block,50,40^-1,0);
    toc
    
    tic
    [mask] = generateRowMask(block,f,phase);
    toc
    
    rowMask = mask < -.85;
    imshow(rowMask,[]);
    out = flattenMaskOverlay(double(bindVec(block)),rowMask,.5,'r');
    figure;
    imshow(out,[]);
    
    
    sig1 = mean(block,1);
    sig1 = zscore(sig1);
    
    sig2 = mean(rowMask,1);
    sig2 = zscore(sig2);
    figure
    plot(sig1,'g')
    hold on
    plot(sig2,'r')
    
    figure;
    imshow(out(:,1:300,:),[]);
    
    
    
    
    
    phun(ridx(r),midx) = mask((end-1)/2,:);
    
    %imshow(phun,[]);
    %drawnow
    r
    numel(ridx)
end
%%
close all
wow = mean(I.*BW,2);
plot(wow.*(wow>.4))
widx = find(wow > .2);
wf = fft(wow(widx));
plot(abs(wf))
bf = 45/size(wf,1);
bm = cos(2*pi*(1:size(I,1))*bf);
bm = repmat(bm',[1 size(I,2)]);
imshow(bm > .95,[]);
out = flattenMaskOverlay(double(I.*BW),bm>.95,.5,'b');
bwf = bm > .95;
%%
phunM = phun < -.95;
%%
out = flattenMaskOverlay(double(bindVec(BW.*I)),logical(BW.*(bwf|phunM)),.7,'r');
%%
close all
plot(mean(subI,1))
%%
close all
phase = angle(sig(:,113));
for r = 1:size(subI,1)
    what(r,:) = cos(2*pi*(1:size(subI,2))*T^-1 + phase(r));
end
imshow(what,[])
msk = what > .8;
out = flattenMaskOverlay(double(subI),msk);
imshow(out,[]);

