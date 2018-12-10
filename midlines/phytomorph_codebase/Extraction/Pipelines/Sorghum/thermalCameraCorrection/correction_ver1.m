I = double(imread('/home/nate/Downloads/DJI_0299_original.tif'));
cI = imread('/home/nate/Downloads/DJI_0299_corected.tif');
nI = adapthisteq(I/max(I(:)),'NumTiles',[11 11]);
nI = nI * max(I(:));
U = mean(nI(:));
U2 = mean(I(:));
nI = nI - U;
nI = nI + U2;
close all
[f1 x1] = ksdensity(cI(:));
[f2 x2] = ksdensity(I(:));
[f3 x3] = ksdensity(nI(:));
%%
close all
figure;
hold on
plot(x1,f1,'r');
plot(x2,f2,'k');
plot(x3,f3,'g');
%%

figure;
imshow([I nI],[]);
%%
close all
figure
imshow(I,[])
%%
figure;
imshow(nI,[])
%%
close all
for e = 3:15
    nI = adapthisteq(I/max(I(:)),'NumTiles',[e e]);
    imshow(nI,[]);
    drawnow
end
%%
s1 = stdfilt(I);
s2 = stdfilt(nI);
s3 = stdfilt(cI);
%%
close all
figure;
imshow(s2,[])
figure;
imshow(s1,[]);
figure;
imshow(s3,[]);
%%
close all
U1 = imfilter(I,fspecial('disk',31),'replicate');
U2 = imfilter(nI,fspecial('disk',31),'replicate');
U3 = imfilter(cI,fspecial('disk',31),'replicate');
figure;
imshow(U1,[]);
figure;
imshow(U2,[]);
figure;
imshow(U3,[]);