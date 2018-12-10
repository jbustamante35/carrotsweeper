%% write out csv data again
FilePath = '/mnt/scratch3/users/nmiller/phytoMorph/Root Morphometrics V2/20-Nov-2013 14:58:27/';
FileList = {};
FileExt = {'csv'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
MT = [];
WT = [];
for i = 1:numel(FileList)
    [pth nm ext] = fileparts(FileList{i});
    pth
    if strcmp(nm,'angle')
        D = csvread(FileList{i});
        D = D(1:300);
        if all(abs(diff(D,1,1)) < 15) & abs(D(1)) < 20
            if ~isempty(strfind(pth,'lip5-1'))
                MT = [MT D];
            elseif ~isempty(strfind(pth,'wt'))
                WT = [WT D];
            end
            plot(D)
            drawnow
            hold on
        end
    end
end
%%
close all
uW = mean(WT,2);
sW = std(WT,1,2).*size(WT,2)^-.5;
uM = mean(MT,2);
sM = std(MT,1,2).*size(MT,2)^-.5;
errorbar(1:numel(uW),uW,sW,'k')
a = axis;
a(1) = 0;
a(2) = numel(uW);
axis(a);
hold on
errorbar(1:numel(uM),uM,sM,'r')
[h p] = ttest2(WT,MT,[],[],[],2);
ax = axes();
plot(1:numel(h),h)
set(ax,'Color','None');
set(ax,'Visible','off');
a = axis;
a(4) = 2;
axis(a)