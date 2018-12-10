inFilePath = '/mnt/snapper/nate/mirror_images/kernelSwelling/Scott/';
oPath = '/mnt/snapper/nate/mirror_images/kernelSwelling/Scott/results/';
mkdir(oPath)
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
cnt = 1;
oData = {};
for e = 1:numel(SET)
    close all
    [swell area para] = generateOutFileBase(SET{e}{1},oPath,1);
    clear Ur Sr Uf Sf p tr LEG A
    for s = 1:12
        try
            % if it can not read then fail because
            D{s} = csvread(strrep(area,'#ROWNUM#',num2str(s)));
            
            % generate for export 
            [pth,nm,ext] = fileparts(area);
            fidx = strfind(nm,'--');
            dateField = strrep(nm(fidx(1)+2:fidx(2)-1),'_','/');
            dateField = [dateField(6:7) '/' dateField(9:10) '/' dateField(3:4)];
            ScannerField = nm(strfind(nm,'Scanner')+7:fidx(3)-1);
            RowField = num2str(s);
            
            % clearn the area data
            [D{s} jAreaData] = cleanSwellData(D{s},2000);
            A{s} = mean(D{s}(1:3,:));
            
            % calc swelling
            D{s} = calcPercentSwelling(D{s}',1);
            % clean swelling data
            [D{s} jSwellData] = cleanSwellData(D{s}',.1);
            D{s} = D{s}';
            % fit the swelling data
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
            %pause(1)
            %waitforbuttonpress
            
            oData{cnt,1} = dateField;
            oData{cnt,2} = ScannerField;
            oData{cnt,3} = RowField;            
            oData{cnt,4} = size(D{s},1);
            
            oData{cnt,5} = mean(p{s}(:,1),1);
            oData{cnt,6} = var(p{s}(:,1),1);
            oData{cnt,7} = mean(p{s}(:,2),1);
            oData{cnt,8} = var(p{s}(:,2),1);            
            oData{cnt,9} = mean(A{s});
            oData{cnt,10} = var(A{s},1);            
            
            for l = 1:size(Uf,2)
                oData{cnt,10+l} = Uf(s,l);
            end
            str = 11 + size(Uf,2);
            
            for l = 1:size(Uf,2)
                oData{cnt,str+l} = Ur(s,l);
            end
            
            cnt = cnt + 1;
            
            
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
oPath = '/mnt/snapper/nate/mirror_images/kernelSwelling/Scott/compiled_results/';
cell2csv([oPath 'Master_results.csv'],oData,',');