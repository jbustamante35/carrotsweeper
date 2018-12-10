file = '/mnt/snapper/Sam Seeds/test3.tif';
file = '/mnt/snapper/Sam Seeds/16-7-11.tif';
file = '/mnt/snapper/Sam Seeds/16-7-12.tif';
I = imread(file);
%%
mut = imcrop(I);
wt = imcrop(I);
%%
imshow(wt,[])
%%
wtOUT = newSeedSizeSimple(wt,100,1);
mutOUT = newSeedSizeSimple(mut,100,0);
%%
close all
[y1 x1] = ksdensity(wtOUT(:,1));
hold on
[y2 x2] = ksdensity(mutOUT(:,1));
plot(x1,y1,'k');
hold on
plot(x2,y2,'r');
[h p] = ttest2(wtOUT,mutOUT,[],[],1)
%%
TH = linspace(-pi,pi,1000)';
wtC = [mean(wtOUT(:,2))/2*cos(TH) mean(wtOUT(:,3))/2*sin(TH)];
mutC = [mean(mutOUT(:,2))/2*cos(TH) mean(mutOUT(:,3))/2*sin(TH)];
close all
plot(wtC(:,1),wtC(:,2),'k');
hold on
plot(mutC(:,1),mutC(:,2),'r');
