inFilePath = '/mnt/snapper/nate/mirror_images/kernelSwelling/Scott/';
oPath = '/mnt/snapper/nate/mirror_images/kernelSwelling/Scott/results/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan for new images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(inFilePath,FileList,FileExt,verbose);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sort SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:numel(SET)
    N = [];
    for img = 1:numel(SET{e})
        [p n ex] = fileparts(SET{e}{img});
        N(img) = str2num(n);
    end
    [N sidx] = sort(N);
    SET{e} = SET{e}(sidx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove those which are in the junk folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:numel(SET)
    fidx = [strfind(SET{e}{1},'junk') strfind(SET{e}{1},'test')];
    
    if ~isempty(fidx)
        rmidx(e) = 1;
    else
        rmidx(e) = 0;
    end
end
SET(find(rmidx)) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find sets with more than 250
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:numel(SET)
    numImages(e) = numel(SET{e});
end
rmidx = numImages < 250;
SET(rmidx) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% browse last image
DS = 10;
PP = [];
MM = [];
MUF = [];
MSF = [];
TA = [];
TR = [];
WD = {};
WDr = {};
for e = 1:numel(SET)
    close all
    [swell area para] = generateOutFileBase(SET{e}{1},oPath,1);
    clear Ur Sr Uf Sf p tr LEG
    for s = 1:12
        try
            D{s} = csvread(strrep(area,'#ROWNUM#',num2str(s)));
            [D{s} jAreaData] = cleanSwellData(D{s},2000);
            A{s} = mean(D{s}(1:3,:));
            D{s} = calcPercentSwelling(D{s}',1);
            [D{s} jSwellData] = cleanSwellData(D{s}',.1);
            D{s} = D{s}';
            [p{s} ep{s} er{s}] = fitSwellBulk(D{s}',size(D{s},2));
            rmidx = er{s} > .8 | isnan(er{s});
            
            % get the junk data
            jD = D{s}(rmidx,:);
            % delete the junk data
            p{s}(rmidx,:) = [];
            ep{s}(rmidx,:) = [];
            er{s}(rmidx) = [];
            D{s}(rmidx,:) = [];
            A{s}(rmidx) = [];

            % average curves
            Ur(s,:) = mean(D{s},1);
            Sr(s,:) = std(D{s},1,1)*size(ep{s},1).^-.5;
            Uf(s,:) = mean(ep{s},1);
            Sf(s,:) = std(ep{s},1,1)*size(ep{s},1).^-.5;        


            % average area
            iA(s) = mean(A{s});

            % average parameters
            P1(s) = mean(p{s}(:,1));
            P2(s) = mean(p{s}(:,2));
            S1(s) = std(p{s}(:,1));
            S2(s) = std(p{s}(:,2));

            % legned
            LEG{s} = num2str(s);

            % trial size
            tr(s) = size(D{s},1);
            
            % raw and fit whole data set
            WD{end+1} = ep{s};
            WDr{end+1} = D{s};
            
            % fly by plots
            plot(D{s}');
            hold on
            plot(ep{s}','k');
            plot(jD','r--');
            hold off
            drawnow
            pause(1)
            %waitforbuttonpress
        catch ME
            ME;
        end
    end
    
    
    
    
    %I = imread(SET{e}{end});
    figure;
    errorbar(Uf(:,1:DS:end)',Sf(:,1:DS:end)');
    legend(LEG);
    
    %figure;
    %imshow(I,[]);    
    %{
    figure;
    bar(RM);
    title('Ratio');
    %}
    figure;
    bar(P1);
    title('Para1');
    
    figure;
    bar(P2);
    title('Para2');
    
    
    % stack across
    PP = [PP ;[P1' P2' S1' S2']];    
    TA = [TA;iA'];
    MUF = [MUF;Uf];
    MSF = [MSF;Sf];
    TR = [TR;tr'];
    
    
    %waitforbuttonpress
end
%%
close all
errorbar(MUF',MSF')
figure;
plot(MUF')
%% look at max
[mv sidx] = max(PP,[],1);
%% 
MAN = [7 5 9 3 8 8 7 7 5 7 6 9 9 4 2 8 6 5 8 4 7 8 7 9 6 4 8 2 8 5 8 8 9 9 8 5 9 9 9 9 9 9 8 6 4 9 9 7 6 1 2 6 1 4];
rmidx = any(isnan(PP),2);
MAN(rmidx) = [];
PP(rmidx,:) = [];
MUF(rmidx,:) = [];
MSF(rmidx,:) = [];
TA(rmidx) = [];
%%
corr(MAN',PP(:,2).*PP(:,1))
corr(MAN',PP)
corr(MAN',TA)
corr(PP,TA)
plot(PP(:,2),MAN,'.')


%%
close all
clear COR PRE

for f = 1:20
    for e = 1:size(MAN,2)iA(rmidx) = [];
        sidx = setdiff(1:numel(MAN),e);
        %[XL,YL,XS,YS,BETA]= plsregress(MM(sidx,:),MAN(sidx)',f);
        %[XL,YL,XS,YS,BETA]= plsregress(MSF(sidx,:),MAN(sidx)',f);
        [XL,YL,XS,YS,BETA]= plsregress([MSF(sidx,:) MUF(sidx,:)],MAN(sidx)',f);
        %[XL,YL,XS,YS,BETA]= plsregress(PP(sidx,:),MAN(sidx)',f);
        PRE(e) = [1 [MSF(e,:) MUF(e,:)]]*BETA;
        %PRE(e) = [1 [PP(e,:)]]*BETA;
    end
    COR(f) = corr(PRE',MAN');
    plot(COR)
    drawnow
    
end
%%
[XL,YL,XS,YS,BETA]= plsregress(MM(sidx,:),MAN(sidx)',f);
