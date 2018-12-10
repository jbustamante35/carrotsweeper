%%
I = bfopen('/home/nate/Downloads/boop/112917_spl02_E4S.czi');
%%
I = bfopen('/home/nate/Downloads/boop/120517/120517_spl01_20x_E4S_phytagel.czi');
%%
I = bfopen('/home/nate/Downloads/boop/120517/120517_spl02_20x_E4S_phytagel.czi');
%%
I = bfopen('/home/nate/Downloads/boop/120117/120117_spl01_20x_E4S.czi');
%%
I = bfopen('/home/nate/Downloads/boop/120117/120117_spl02_20x_E4S.czi');
%%
I = bfopen('/home/nate/Downloads/boop/ecoliData/080617/080617_spl07_E1S_tower of the dead_20x.czi');
%% look at single cell traces - ecoli
M1 = [];
M2 = [];
% load the gfp and grayscale images
for e = 1:2:size(I{1},1)
    M1 = cat(3,M1,I{1}{e,1});
    M2 = cat(3,M2,I{1}{e+1,1});
    e
end
M1 = double(M1);
%%
for e = 1:size(M1,3)
    M1(:,:,e) = imfilter(M1(:,:,e),fspecial('average',[3 3]),'replicate');
end
%%
sz1 = size(M1);
S1 = reshape(M1,[prod(sz1(1:2)) sz1(3)]);
%%
% smooth the data
sS1 = imfilter(S1,fspecial('average',[1 11]),'replicate');
BK = imfilter(S1,fspecial('average',[1 101]),'replicate');
sS1 = sS1 - BK;
%%
close all
plot(sS1(5000,:))
hold on
plot(S1(5000,:))
%% look at smoothed movie
close all
J = reshape(S1,sz1);
for e = 1:size(J,3)
    imshow(J(:,:,e),[])
    drawnow
end
%%
close all
[c r V] = impixel(M1(:,:,1)/255);
sig = J(r,c,:);
plot(squeeze(sig))
%% look at many
close all
%plot(sS1(1:1000:end,:)');
ft = abs(fft(sS1-mean(sS1,2),[],2));
uft = mean(ft,1);
plot(uft);

sig = ft(:,1:20);

NC = 5;
kidx = kmeans(sig,NC);
R = [];
U = [];
for u = 1:NC
    fidx = find((kidx==u));
    
    
    R(u,:) = sS1(fidx(300),:);
    
    U(u,:) = mean(sS1(fidx,:));
    
    plot(mean(ft(fidx,:),1))
    hold on
end
%%
close all
figure;
plot(U')
title('mean');
figure;
plot(R');
%%
close all
K = reshape(kidx,[sz1(1:2)]);
RGB = label2rgb(K);
imshow(RGB,[]);
%%
for e = 1:NC
    close all
    figure
    imshow(K==e,[]);
    R = regionprops(K==e,'Area','PixelIdxList');
    VEC = [];
    for r = 1:numel(R)
        VEC = [VEC;mean(sS1(R(r).PixelIdxList,:),1)];
    end
    
    figure;
    plot(VEC(1:10,:)');
    
    figure
    plot(U(e,:))
    waitforbuttonpress
end
%%
imageSize = [528 512];
%% load whole movie
M1 = [];
M2 = [];
% load the gfp and grayscale images
for e = 1:2:size(I{1},1)
    M1 = cat(3,M1,I{1}{e,1});
    M2 = cat(3,M2,I{1}{e+1,1});
    e
end
% convert to double
M1 = double(M1);
M2 = double(M2);
sz = size(M1);
% reshape data
R1 = reshape(M1,[prod(sz(1:2)) sz(3)]);
R2 = reshape(M2,[prod(sz(1:2)) sz(3)]);

% get the pixel level baseline
in1 = imfilter(R1,fspecial('average',[1 101]),'replicate');
in2 = imfilter(R2,fspecial('average',[1 101]),'replicate');
% smooth the data
R1 = imfilter(R1,fspecial('average',[1 11]),'replicate');
R2 = imfilter(R2,fspecial('average',[1 11]),'replicate');close all

% subtract pixel level baseline from smoothed data
my1 = R1 - in1;
my2 = R2 - in2;

% get the baseline over frames of the smoothed curve
U1 = mean(R1,1);
U2 = mean(R2,1);

% normalize the b
N1 = bindVec(U1);
N2 = bindVec(U2);

% extract the baseline
baseline1 = imfilter(U1,fspecial('average',[1 101]),'replicate');
baseline2 = imfilter(U2,fspecial('average',[1 101]),'replicate');
D1 = bsxfun(@minus,double(R1),baseline1);
D2 = bsxfun(@minus,double(R2),baseline2);
F1 = mean(D1,1);
F2 = mean(D2,1);
plot(F1);hold on
plot(F2);
%% plot amp
close all
plot(U1,'g')
hold on
figure;
plot(U2,'k')
%% view "raw" data
close all
p1 = min(my1(:));
p2 = max(my1(:));

p1 = std(my1(:));


p1 = std(my1(:));
p2 = std(my2(:));

h1 = figure;
h2 = figure;
LK1 = [];
LK2 = [];
for e = 1:size(my1,2)
    W1 = reshape(my1(:,e),imageSize);
    W2 = reshape(my2(:,e),imageSize);
    W1 = imfilter(W1,fspecial('average',[21 21]),'replicate');
    W2 = imfilter(W2,fspecial('average',[21 21]),'replicate');
    l1(e) = mean(W1(:));
    l2(e) = mean(W2(:));
    
    
    fidx = find(MSK);
    LK1(e) = mean(W1(fidx));
    LK2(e) = mean(W2(:));
    W2 = W2 * p2^-1;
    W2 = W2 * p1;
    WT = cat(3,W1,zeros(size(W1)),W2);
    figure(h1);
    imshow(WT,[-p1 p2]);
    figure(h2)
    plot(LK1,'g');
    hold on
    plot(LK2,'k');
    hold off
    drawnow
   
    e
end
%% find areas of occilations
close all
ST1 = std(my1,[],2);
ST1 = reshape(ST1,imageSize);
MSK = ST1 > graythresh(ST1);
imshow(MSK,[]);
%%
close all
for e = 1:20:size(M2,3)
    tmp = M2(:,:,e);
    tmp = imfilter(tmp,fspecial('average',[5 5]),'replicate');
    contour(tmp)
    %imshow(tmp,[50 160]);
    drawnow
end
%%
close all
J1 = abs(fft(bsxfun(@minus,D1,mean(D1,2)),[],2));
A1 = angle(fft(bsxfun(@minus,D1,mean(D1,2)),[],2));
plot(mean(J1,1),'g')
freq = 9;
C1 = reshape(J1(:,freq),imageSize);
hold on
J2 = abs(fft(bsxfun(@minus,D2,mean(D2,2)),[],2));
A2 = angle(fft(bsxfun(@minus,D2,mean(D2,2)),[],2));
plot(mean(J2,1),'k')
C2 = reshape(J2(:,freq),imageSize);

A1 = A1(:,freq);

A2 = A2(:,freq);



figure;
imshow(C1,[]);
figure;
imshow(C2,[]);



sC1 = imfilter(C1,fspecial('average',[21 21]),'replicate');
sC2 = imfilter(C2,fspecial('average',[21 21]),'replicate');
figure;
imshow(sC1,[]);
figure;
imshow(sC2,[]);






%%
figure
imshow(cat(3,bindVec(sC1),zeros(size(sC1)),bindVec(sC2)),[]);
%%
close all
S1 = [];
S2 = [];
h1 = figure;
h2 = figure;
for t = 1:size(D1,2)
    %{
    T1 = C1(:).*sin(t*2*pi/93 + A1);
    T2 = C2(:).*sin(t*2*pi/93 + A2);
    T1 = reshape(T1,imageSize);
    T2 = reshape(T2,imageSize);
    %}
    d1 = D1(:,t);
    %d1 = d1 - mean(d1);
    d2 = D2(:,t);
    %d2 = d2 - mean(d2);
    d1 = reshape(d1,imageSize);
    d2 = reshape(d2,imageSize);
    
    
    
    
    %{
    d1 = d1 - min(D1(:));
    d2 = d2 - min(D2(:));
    d1 = d1 * max(D1(:))^-1;
    d2 = d2 * max(D2(:))^-1;
    %}
    figure(h1)
    imshow(cat(3,d1,zeros(size(d1)),d2),[]);
    
    %{
    T1 = T1 .* max(abs(C1(:)))^-1;
    T2 = T2 .* max(abs(C2(:)))^-1;
    figure(h2);
    imshow(cat(3,T1,zeros(size(T1)),T2),[]);
    %}
    S1(t) = mean(T1(:));
    S2(t) = mean(T2(:));
    t
    
    
    
    s1(t) = mean(d1(:));
    s2(t) = mean(d2(:));
    drawnow
end
%%
close all
for e = 1:size(D1,1)
    vec1 = D1(e,:);
    vec2 = D2(e,:);
    vec1 = imfilter(vec1,fspecial('average',[1 51]),'replicate');
    vec2 = imfilter(vec2,fspecial('average',[1 51]),'replicate');
    plot(vec1 - mean(vec1),'g');
    hold on
    plot(vec2-mean(vec2),'k');
    hold off 
    axis([1 numel(vec1) -2 2]);
    drawnow
end



%%
angle = pi/2;
angle = linspace(0,pi,1);
A = [];
 close all
h1 = figure;
h2 = figure;
for a = 1:numel(angle)
    [R T] = ndgrid(linspace(-30,30,60),linspace(angle(a),angle(a),1));
    X1 = R.*sin(T);
    X2 = R.*cos(T);
    X1 = X1 + r;
    X2 = X2 + c;
    %
   
    W = [];
    disp = 0;
    cnt = 1;
    G = [];
    for e = 1:2:size(I{1},1)
        tmpI = I{1}{e,1};
        tmpI = imfilter(double(tmpI)/255,fspecial('disk',11),'replicate');
        tmpI = tmpI - offset(cnt);
        vec(cnt) = mean(tmpI(:));
        figure(h1);
        contour(tmpI,linspace(0,.1,20));
        axis off
        drawnow
        G(:,:,:,cnt) = frame2im(getframe());
        figure(h2);
        plot(vec);
        drawnow;
        cnt = cnt + 1;
        what = ba_interp2(tmpI,X2,X1);
        W = cat(2,W,what);
        if disp
            if e == 1
                imshow(tmpI,[]);
                hold on
                plot(X2,X1,'.')
                hold off
                drawnow
                e
            end
        end
    end
    A = cat(3,A,W);
    a
end
%%
cnt = 1;
for e = 1:2:size(I{1},1)
    tmpI = I{1}{e,1};
    tmpI = imfilter(double(tmpI)/255,fspecial('disk',11),'replicate');
    vec(cnt) = mean(tmpI(:));
    cnt = cnt + 1;
end

close all
offset = imfilter(vec,fspecial('average',[1 131]),'replicate');
plot(vec-offset)
%%
close all
figure
for e = 1:size(A,3)
    mesh(A(:,:,e));
    view([0 90]);
    title(num2str(angle(e)*180/pi));
    drawnow
end
%%
close all
h1 = figure;
h2 = figure;
svec = imfilter(vec,fspecial('average',[1 41]),'replicate');
for e = 1:size(G,4)
    figure(h1)
    imshow(G(:,:,:,e)/255,[]);
    figure(h2);
    plot(svec(1:e))
    
    drawnow
end
%%
%%

svec = imfilter(vec-offset,fspecial('average',[1 61]),'replicate');
close all
plot(svec/max(svec(:)))
hold on
upper = imdilate(svec,ones(1,51)) == svec;
lower = imerode(svec,ones(1,51)) == svec;
plot(upper,'r')
plot(lower,'k')
uidx = find(upper);
lidx = find(lower);
waitforbuttonpress
close all
level = linspace(-.02,.06,5);
for e = 1:numel(uidx)
    
    
    UtmpI = I{1}{2*(uidx(e)-1)+1,1};
    UtmpI = imfilter(double(UtmpI)/255,fspecial('disk',11),'replicate');
    UtmpI = UtmpI - offset(uidx(e));
    
    
    
    
    LtmpI = I{1}{2*(lidx(e)-1)+1,1};
    LtmpI = imfilter(double(LtmpI)/255,fspecial('disk',11),'replicate');
    LtmpI = LtmpI - offset(lidx(e));
    
    
    contour(UtmpI,level);
    hold on
    contour(LtmpI,level);
    waitforbuttonpress
    hold off
    
    
   % drawnow
    %pause(.5)
    
end
%%
close all
[J,BOX] = imcrop(I{1}{1,1},[]);
close all
for e = 1:8:size(I{1},1)
    tmpI = I{1}{e,1};
    tmpI = imcrop(tmpI,BOX);
    tmpI = imfilter(tmpI,fspecial('average',[41 41]),'replicate');
    imshow(tmpI,[]);
    drawnow
end