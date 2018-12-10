fn = '/home/nate/Downloads/170117 LC1-107-2.xlsx';
[n,~,d] = xlsread(fn,2);
n(:,1:2) = [];
n(1:2,:) = [];
%% dig for data
FilePath = '/home/nate/Downloads/Alex/';
FileList = {};
FileExt = {'xlsx'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% read sheets
clear d
parfor e = 1:numel(FileList)
    [p,nm,ext] = fileparts(FileList{e});
    nm = strrep(nm,' LC1-','');
    fidx = strfind(nm,'-');
    rep = nm(fidx(1)+1);
    nm = nm(1:fidx(1)-1);
    
   
    
    [n,~,~] = xlsread(FileList{e},2);
    n(:,1:2) = [];
    n(1:2,:) = [];
    wellName = 1:size(n,2);
    chemName = [];
    for i = 1:numel(wellName)
        chemName(i) = str2num([num2str(wellName(i)) nm]);
    end
    d(e).data = n;
    d(e).d1 = n(1:13,:);
    d(e).d2 = n(13:194,:);
    d(e).d3 = n(195:end,:);
    d(e).name = chemName
    d(e).rep = ones(1,size(n,2))*rep;
    
    e
end
%% stack all (2)
clear n
n = [d.data];
[nn] = normalizeData(n);
close all
rmidx = any(diff(nn,1,1)  == 0,1);
chemName = [d.name];
chemRep = [d.rep];
%% water index (4)
close all
widx = mod(0:(size(n,2)-1),12);
widx = 1:13:size(n,2);
%% plot water
close all
plot(n(:,widx))
%% decompose all data (4.5)
[S C U E L ERR LAM] = PCA_FIT_FULL_T(nn,5);
%% get and normalize water data (5)
wD = nn(:,widx);
wC = C(:,widx);
[nwD] = normalizeData(wD);
rwi = rmidx(widx);
wD(:,rwi) = [];
%% get and normalize non water data (6)
midx = setdiff(1:size(nn,2),widx);
mC = C(:,midx);

chemName = chemName(midx);
chemRep = chemRep(midx);
mD = nn(:,midx);
rmi = rmidx(midx);

mD(:,rmi) = [];
chemName(rmi) = [];
repName(rmi) = [];


[nmD] = normalizeData(mD);
d1 = nmD(1:13,:);
d2 = nmD(13:194,:);
d3 = nmD(195:end,:);
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL_T(d1,5);
[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL_T(d2,5);
[S3 C3 U3 E3 L3 ERR3 LAM3] = PCA_FIT_FULL_T(d3,5);

%% kmeans
kidx = kmeans(mC',3);
for k = 1:5
    sub{k} = nmD(:,kidx==k);
end
%% keams 2
close all
kidx = kmeans(C2',3);
grp = kidx;
CL = {'r' 'b' 'g' 'c' };
for k = 1:3
    % pull out the kth group from the first cluster and cluster
    tidx = kmeans(C3(:,kidx==k)',3);
    fidx = find(kidx==k);
    u0 = mean(d1(:,kidx==k),2);
    s0 = std(d1(:,kidx==k),1,2)*sum(kidx==k)^-.5;
    u1 = mean(d2(:,kidx==k),2);
    s1 = std(d2(:,kidx==k),1,2)*sum(kidx==k)^-.5;
    
    
    
    for u = 1:3
        % find the 
        grp(fidx(tidx==u),2) = u;
        u2 = mean(d3(:,fidx(tidx==u)),2);
        s2 = std(d3(:,fidx(tidx==u)),1,2)*sum(tidx==u)^-.5;
        %plot([u0;u1;u2]);
        %errorbar([u0;u1;u2],[s0;s1;s2]);
        grpU{k,u} = [u0;u1;u2];
        grpS{k,u} = [s0;s1;s2];
        plot([u0;u1;u2],CL{u});
        hold on
        plot([u0;u1;u2]+10*[s0;s1;s2],[CL{u} '--']);
        plot([u0;u1;u2]-10*[s0;s1;s2],[CL{u} '--']);
        LEG{u} = num2str(sum(tidx==u));
    end
    
    
    
    se = std(wD,1,2).*size(wD,2).^-.5;
    plot(mean(wD,2),'k');
    plot(mean(wD,2)+10*se,'k--');
    plot(mean(wD,2)-10*se,'k--');
    
    %errorbar(mean(wD,2),,'k');
    legend(LEG);
    waitforbuttonpress
    hold off
end
%% find those that are in same group
close all

UQ = unique(chemName);
for u = 1:numel(UQ)
    fidx = find(chemName == UQ(u));
    tG = grp(fidx,:);
    if size(tG,1) == 2
        SIM(fidx) = all(tG(1,:) == tG(2,:));
        whichG(fidx) = (tG(1,1)-1)*3 + (tG(2,1)-1) + 1;
        Y(fidx) = UQ(u);
        subG(fidx,:) = tG;
    else
        SIM(fidx) = 0;
        whichG(fidx) = 0;
        Y(fidx) = UQ(u);
    end
end
FUN = whichG.*SIM;
%% which are the groups
close all
sidx = find(SIM);
g = grp(sidx,:);
UQg = unique(g,'rows');
for u = 1:numel(UQg)
    plot(grpU{UQg(u,1),UQg(u,2)});
    
    hold on
    
end
%% plot water vs non water
close all
plot(mean(wD,2),'b');
hold all
plot(mean(mD,2),'r');
for e = 1:numel(sub)
    plot(mean(sub{e},2),'k')
end
%% look for water pattern
close all
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL_T(nwD,3);
plot(wS);
figure
plot3(wC(1,:),wC(2,:),wC(3,:),'.');
figure;
plot(wU)
%% split into water and non water

%% classify data
close all
[aS aC aU aE aL aERR aLAM] = PCA_FIT_FULL_T(nn,3);
plot(aS);
figure
plot3(aC(1,:),aC(2,:),aC(3,:),'.');
figure;
plot(aU)
%% 
close all
plot(n);
%%
close all
[S C U E L ERR LAM] = PCA_FIT_FULL_T(n,6);
plot(S);
%% split too
d1 = n(1:13,:);
d2 = n(13:194,:);
d3 = n(195:end,:);
%% look at raw and fit
close all
for e = 1:size(S,2)
    plot(n(:,e),'k');
    hold on
    plot(S(:,e),'r');
    hold off 
    waitforbuttonpress
end
%% break down d1
close all
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL_T(d1,1);
plot(S1);
%% break down d2
close all
[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL_T(d2,3);
plot(S2);
for e = 1:size(S2,2)
    plot(d2(:,e),'k');
    hold on
    plot(S2(:,e),'r');
    hold off
    waitforbuttonpress
end
%% break down d3
close all
[S3 C3 U3 E3 L3 ERR3 LAM3] = PCA_FIT_FULL_T(d3,3);
plot(S3);
for e = 1:size(S3,2)
    plot(d3(:,e),'k');
    hold on
    plot(S3(:,e),'r');
    hold off
    waitforbuttonpress
end
%%
