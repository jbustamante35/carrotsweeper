% read data for emergnece
clear d
d1 = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/20170131_Camera3.csv');
d2 = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/20170131_Camera4.csv');

d(1).experimentName = '20170131_Camera3';
d(2).experimentName = '20170131_Camera4';
d(1).data = cell2mat(d1);
d(2).data = cell2mat(d2);
%% emergnece master list
ml = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/EmergenceAssay_Master_list.csv');
%% auto generate d structe
FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/return/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(FileList)
    fileName = FileList{e};
    fidx = strfind(fileName,'_');
    [pth,fn,ext] = fileparts(FileList{e});
    tmp = readtext(FileList{e});
    
    d(e).data = cell2mat(tmp);
    d(e).experimentName = fn;
end
%% other data

em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/NAM_parents_ear_map_xyCoords(1).csv');
%% for swelling

FilePath = '/mnt/snapper/nate/mirror_images/maizeData/jgustin/kernel_swellData/return/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
kp = [];
pkp = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e},'swell--'))
        kp = [kp e];
    end
     if ~isempty(strfind(FileList{e},'para--'))
        pkp = [pkp e];
    end
end
pFileList = FileList(pkp);
FileList = FileList(kp);

close all
im = readtext('/home/nate/Downloads/imbibition_assay_master_list(1).csv');
LLI = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K'};
import java.util.Map.*;
import java.util.HashMap;
SM = HashMap();
for e = 2:size(im,1)
    earKey = im{e,6}
    scannerKey = im{e,4};
    arrayName = im{e,10};
    arrayName = ['; ' arrayName ';'];
    fidx = strfind(arrayName,';');
    fidx = strfind(arrayName,';');
    for l = 1:(numel(fidx)-1)
        snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
        sidx = strfind(snip,':');
        letterRange1 = snip(2);
        letterRange2 = snip(sidx(1)+1);
        numberRange1 = str2num(snip(3:(sidx(1)-1)));
        numberRange2 = str2num(snip((sidx(1)+2):end));
        
        rowN = find(strcmp(LLI,letterRange1));
        for l = 1:numel(FileList)
            if ~isempty(strfind(FileList{l},['--' num2str(rowN)])) & ~isempty(strfind(FileList{l},scannerKey))
                sd = csvread(FileList{l});
                pd = csvread(pFileList{l});
                [keep] = manfredFilter(sd');
                plot(sd,'k')
                hold on
                plot(sd(:,keep),'r')
                hold off
                drawnow
                %waitforbuttonpress
                earMean = mean(sd(:,keep),2);
                earMean = mean(pd(keep,2));
                SM.put(earKey,earMean);
            end
        end
    end
end
%% parse ear map
cameraName = 8;
posML2 = 9;
import java.util.Map.*;
import java.util.HashMap;
EM = HashMap();
for e = 2:size(em,1)
    if ~strcmp(em{e,cameraName},'NA')
        key = ['20170131_' lower(em{e,cameraName}) '-' lower(em{e,posML2}(1)) '-' lower(em{e,posML2}(2:end))];
        value = em{e,12} - em{e,14};
        LEN = em{e,15} - em{e,14};
        value = value/LEN;
        if strcmp(key,'20170131_camera4-d-g12')
            value
        end
        if isempty(value)
            key;
        end
        key;
        EM.put(key,value);
    end
end
%% make labels
import java.util.Map.*;
import java.util.HashMap;
LT = HashMap();

LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
HL = {'d' 'p'};
LABELS = {};
for e1 = 1:numel(HL)
    for e2 = 1:numel(LL)
        for e3 = 1:numel(NL)
            LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
        end
    end
end
for e = 1:numel(d)
    for f = 1:size(d(e).data,2)
        key = lower([d(e).experimentName '-' LABELS{f}]);
        value = d(e).data(:,f);
        LT.put(key,value);
    end
end
%% align data from mastersheet
iplantName = 27+1;
hcName = 20+1;
posName = 10+1;
genoName = 7+1;
import java.util.Map.*;
import java.util.HashMap;
GT = HashMap();
MT = HashMap();
GV = HashMap();
SM2 = HashMap();
for e = 2:size(ml,1)
    for f = 1:numel(d)
        if ~isempty(strfind(ml{e,iplantName},d(f).experimentName))
            
            arrayName = ml{e,posName};
            arrayName = ['; ' arrayName ';'];
            fidx = strfind(arrayName,';');
            for l = 1:(numel(fidx)-1)
                snip = arrayName((fidx(l)+1):(fidx(l+1)-1));
                sidx = strfind(snip,':');
                letterRange1 = snip(2);
                letterRange2 = snip(sidx(1)+1);
                numberRange1 = str2num(snip(3:(sidx(1)-1)));
                numberRange2 = str2num(snip((sidx(1)+2):end));
                gidx = strfind(ml{e,iplantName},filesep);
                cn = numberRange1:numberRange2;
                culN = ml{e,6};
                for rn = 1:numel(LL)
                    for c = 1:numel(cn)
                        key = [ml{e,iplantName}((gidx(end)+1):end) '-' lower(ml{e,hcName}(1)) '-' LL{rn} num2str(cn(c))];
                        flag  = ~strcmp(ml{e,6},'15A-6197-02');
                        
                        key = lower(key);
                        value = LT.get(key);
                        deltaT = find(value);
                        genoTmp = ml{e,genoName};
                        %posTmp = EM.get(key);
                        vec = GT.get(genoTmp);
                        vec = [vec value(1:250)];
                        if flag
                            GT.put(genoTmp,vec);
                        end
                        
                        
                        if isempty(deltaT)
                            deltaT = inf;
                        end
                        %corrV = [posTmp,deltaT(1)];
                        %if flag
                        %    MT.put(key,corrV);
                        %end
                        
                        vec = GV.get(genoTmp);
                        %vec = [vec corrV'];
                        %if flag
                        %    GV.put(genoTmp,vec);
                        %end
                        
                        vec = SM2.get(culN);
                        vec = [vec;deltaT(1)];
                        if flag
                            SM2.put(culN,vec);
                        end
                        
                    end
                end
            end
            lower(ml{e,hcName}(1));
        end
    end
end
%%
close all
keys = MT.keySet;
itr = keys.iterator();
STACK = [];

while itr.hasNext()
    key = itr.next();
    value = MT.get(key);
    
    if ~isinf(value(2))
        plot(value(1),value(2),'.')
        STACK = [STACK;value'];
        hold on
    end
end
%%
close all
keys = SM.keySet;
itr = keys.iterator();
while itr.hasNext()
    key = itr.next();
    val1 = SM.get(key);
    val2 = SM2.get(key);
    val2(isinf(val2)) = [];
    plot(mean(val2),val1,'.')
    hold on
end
%%
close all
keys = GT.keySet;
itr = keys.iterator();
%h1 = figure;
%h2 = figure;
cnt = 1;
mPDF = [];
iPDF = [];
while itr.hasNext()
    %for e = 1
    %key = itr.next();
    key = itr.next();
    vec = GT.get(key);
    %{
    %whole = GV.get(key);
    
    %nonGerm = isinf(whole(2,:));
    %figure(h2);
    %sum(nonGerm)
    f = ksdensity(whole(1,nonGerm),linspace(0,1,1000));
    f = f / sum(f);
    plot(linspace(0,1,1000),f)
    hold all
    plot(whole(1,:),ones(1,numel(nonGerm))*cnt*.001,'.');
    %}
    
    
    
    
    
    %{
    figure(h1);
    
    
    
    
    whole(:,isinf(whole(2,:))) = [];
    plot(whole(1,:),whole(2,:),'.');
    hold all
    %}
    
    no_germN = sum(all(vec==0,1));
    germN = size(vec,2) - no_germN;
    germP(cnt,:) = [germN no_germN]/size(vec,2);
    
    ipdf = gradient(imfilter(mean(vec,2),fspecial('average',[5 1]),'replicate'));
    pdf = imfilter(mean(vec,2),fspecial('average',[5 1]),'replicate')
    ipdf = ipdf / 30  * 60 ;
    ipdf(ipdf < 0) = 0;
    xlab = 1:numel(ipdf);
    xlab = (xlab*30/60 + 80)/24;
    
    
    [J xval] = min(abs(pdf - mean(pdf)));
    [para{cnt}] = fminsearch(@(X)mySigmoid_ver0(xlab',X,pdf),[pdf(end) .5 xlab(xval)]); 
    [e,yp] = mySigmoid_ver0(xlab',para{cnt});
    ypg = gradient(yp);
    %{
    type = 'logit';
    type = 'probit';
    n = size(vec,2)*ones(size(pdf));
    xlab = 1:size(vec,1);
    toFit = sum(vec,2);
    fidx = find(toFit);
    [logitCoef,dev] = glmfit(xlab(fidx(1):end)',[toFit(fidx(1):end) n(fidx(1):end)],'binomial','link',type);
    logitFit = glmval(logitCoef,xlab',type,'size', n);
    plot(xlab,logitFit);
    %}
    
    
    figure;
    [AX,H1,H2] = plotyy(xlab,ipdf*100,xlab,pdf);
    hold(AX(2))
    hold(AX(1))
    plot(AX(2),xlab',yp,'r')
    plot(AX(1),xlab',ypg*100/ 30  * 60 ,'r')
    set(AX(2),'YTick',linspace(0,1,11));
    axis(AX(2),[3 9 0 1]);
    axis(AX(1),[3 9 0 10]);
    set(AX(1),'YTick',linspace(0,10,11));
    set(get(AX(1),'Ylabel'),'String','percent germ per hour') 
    set(get(AX(2),'Ylabel'),'String','total percent germ')
    xlabel('days after planting');
    
    LEG{cnt} = key;
    mPDF(:,cnt) = yp;
    iPDF(:,cnt) = ypg;
    %title([key '-' num2str(size(vec,2))])
    title([key])
    %{
     hold on
    figure;
    plot(xlab',yp)
    hold on
    plot(xlab',pdf,'r')    
    %}
    cnt = cnt + 1;
end
%%
CL = {'r' 'g' 'b' 'm' 'c' 'k' 'r--' 'b--' 'g--'};
h = figure;
ax(1) = axes();
ax(2) = axes('YAxisLocation','right');
for e = 1:size(mPDF,2)
    plot(ax(1),xlab',mPDF(:,e),CL{e});
    hold(ax(1),'on');
    plot(ax(2),xlab',100*iPDF(:,e)/ 30  * 60 ,CL{e});
    hold(ax(2),'on');
end

axis(ax(1),[3 9 0 1]);
set(ax(1),'YTick',linspace(0,1,11));
axis(ax(2),[3 9 0 10]);
set(ax(2),'YTick',linspace(0,10,11));
set(ax(2),'YAxisLocation','right');
set(ax(2),'Color','none');
ylabel(ax(1),'Percent Germination')
ylabel(ax(2),'Percent/hour Germination')
xlabel(ax(1),'Time (days)')
legend(LEG)
%%

figure;
bar(germP,'stack')



