%% load ALL POP specdata from database
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maizelj2','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maizelj2');
q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...                
                'JOIN predictions ' ...
                'on kernels.id = predictions.kernel_id ' ...                
                'JOIN averageweightspectra_vw ' ...
                'ON kernels.id=averageweightspectra_vw.kernel_id ' ...                
                'JOIN files ' ...
                'on kernels.img_gray_front_fid = files.id ' ...
                'WHERE population_lines.type = ' '''NC-350_RILs'' OR population_lines.type = ''NAM_parents'' OR population_lines.type = ''composition_mutants'' OR population_lines.type = ''seedling_phenotyping_widiv'' '];

q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...                              
                'JOIN averageweightspectra_vw ' ...
                'ON kernels.id=averageweightspectra_vw.kernel_id ' ...
                'JOIN files ' ...
                'on kernels.img_gray_front_fid = files.id ' ...
                'WHERE population_lines.type = ' '''NC-350_RILs'' OR population_lines.type = ''NAM_parents'' OR population_lines.type = ''composition_mutants'' OR population_lines.type = ''seedling_phenotyping_widiv'' '];

cursor = exec(conn, q);
cursor = fetch(cursor);
fieldString = columnnames(cursor,1);
results = cursor.Data;
%{
img_bw_front_fid    | integer          | 
 img_bw_side_fid     | integer          | 
 img_bw_top_fid      | integer          | 
 img_gray_front_fid  | integer          | 
 img_gray_side_fid   | integer          | 
 img_gray_top_fid    | integer          | 
 img_color_front_fid
%}
%% SPLIT ALL
clear A
fidx = find(strcmp(fieldString,'plate_name'));
A.plateName = results(:,fidx);
fidx = find(strcmp(fieldString,'plate_position'));
A.position = results(:,fidx);
fidx = find(strcmp(fieldString,'kernel_id'));
A.kernel_id = results(:,fidx(1));
fidx = [find(strcmp(fieldString,'wl_907')) find(strcmp(fieldString,'wl_1688'))];
A.specData = cell2mat(results(:,fidx(1):fidx(2)));
%A.prediction = cell2mat(results(:,35:40));
A.genoType = results(:,16);
A.POP = {};
for e = 1:size(results,1)
    A.POP{end+1} = char(results{e,3});
end
A.imageTOP = results(:,end);
[kernelVec genoVec popVec kernelID rkernelID] = translateWellNames(A);
%save('/mnt/scratch1/phytoM/flashProjects/maize/specKey2.mat','kernelVec','genoVec','popVec','kernelID');
%% NEW-NEW FRONT 
%ext = '_S2.tiff';
ext = '_F.tiff';
%ext = '_T.tiff';
close all
%Pnew = zeros(numel(f.imageTOP),2);
imageList = A.imageTOP;

clear tmp;
MASTERI = {};
MASTERBOX = {};
MASTERBOXl = {};
qBOX = {};
disp = 0;
%widx = find(strcmp(A.genoType,'W22^ACR'));
%imageList = imageList(widx);
%imageList = imageList(randperm(numel(imageList)));
parfor e = 1:numel(imageList)
    try 
        tmp = readKernelImageFile(imageList{e},ext);
        if ~isempty(tmp)
            tmpf = stdfilt(log(double(tmp)),ones(5));
            tmpf = imfilter(tmpf,ones(11),'replicate');
            [d1 d2] = gradient(tmpf);
            d1 = abs(d1);
            s2 = std(double(d1),1,2);
            s2 = log(s2);
            s2 = imfilter(s2,ones(11,1),'replicate');
            s2 = bindVec(s2);
            R = regionprops(s2 > .4,'Area','PixelIdxList','Centroid');
            cen = [R.Centroid];
            cen = cen(2:2:end);
            [~,sidx] = min(abs(cen - size(tmp,1)/2));
            s2m = zeros(size(s2));
            %[~,sidx] = sort([R.Area]);

            s2m(R(sidx(end)).PixelIdxList) = 1;
            BOX = [0 min(find(s2m)) size(tmp,2) max(find(s2m))-min(find(s2m))];
            BOX = [0 0 size(tmp,2) max(find(s2m))];
            lBOX = BOX;
            lBOX(4) = lBOX(4) - 50;

            MASTERI{e} = tmp;
            MASTERBOX{e} = BOX;
            MASTERBOXl{e} = lBOX;


            ltmp = imcrop(tmp,lBOX);
            ltmp = stdfilt(ltmp,ones(11));
            ltmp = imfilter(ltmp,fspecial('average',5),'replicate');
            s1 = std(ltmp,1,1);
            s1 = bindVec(log(s1+.1));

            R = regionprops(s1 > .1,'Area','PixelIdxList');



            s1 = zeros(size(s1));
            [~,sidx] = sort([R.Area]);
            s1(R(sidx(end)).PixelIdxList) = 1;
            wididx = find(s1);
            ibk(e) = backGroundKernel(tmp);
            mBOX = [min(wididx) 0 sum(s1) lBOX(4)];
            MASTERBOXcrop{e} = mBOX;
            jtmp = imcrop(ltmp,mBOX);
            s2 = std(jtmp,1,2);
            upper = find(s2 > 1);
            upper = min(upper);
            MBOX = mBOX;
            MBOX(2) = upper;
            MBOX(4) = BOX(4) - upper;
            qBOX{e} = MBOX;
            kid(e) = A.kernel_id{e};
            POPU{e} = A.POP{e};
         %{
            if disp
                imshow(tmp,[])
                title(imageList{(e)});
                hold on

                rectangle('Position',MASTERBOX{e},'EdgeColor','r')
                rectangle('Position',MASTERBOXl{e},'EdgeColor','g')
                rectangle('Position',MASTERBOXcrop{e},'EdgeColor','y')
                rectangle('Position',qBOX{e},'EdgeColor','c')

                hold off
                drawnow
                pause(.3)
            end
            %}
            


        else
            e;
            numel(imageList);
            hello = 1;
        end
     
    e
    numel(imageList)
        
        
        
        
        
    catch ME2
        %MASTERI{e} = [];
        %MASTERBOX{e} = [];
    
    end
    
end
%% NEW-NEW SIDE
ext = '_S2.tiff';
%ext = '_F.tiff';
close all
%Pnew = zeros(numel(f.imageTOP),2);
imageList = A.imageTOP;

clear tmp;
sMASTERI = {};
sMASTERBOX = {};
sMASTERBOXl = {};
sqBOX = {};
disp = 1;
%widx = find(strcmp(A.genoType,'W22^ACR'));
%imageList = imageList(widx);
%imageList = imageList(randperm(numel(imageList)));
parfor e = 1:numel(imageList)
    e
    
    try
        tmp = readKernelImageFile(imageList{e},ext);
        if ~isempty(tmp)
            tmpf = stdfilt(log(double(tmp)),ones(5));
            tmpf = imfilter(tmpf,ones(11),'replicate');
            [d1 d2] = gradient(tmpf);
            d1 = abs(d1);
            s2 = std(double(d1),1,2);
            s2 = log(s2);
            s2 = imfilter(s2,ones(11,1),'replicate');
            s2 = bindVec(s2);
            R = regionprops(s2 > .43,'Area','PixelIdxList','Centroid');
            cen = [R.Centroid];
            cen = cen(2:2:end);
            [~,sidx] = min(abs(cen - size(tmp,1)/2));
            s2m = zeros(size(s2));
            %[~,sidx] = sort([R.Area]);

            s2m(R(sidx(end)).PixelIdxList) = 1;
            BOX = [0 min(find(s2m)) size(tmp,2) max(find(s2m))-min(find(s2m))];
            BOX = [0 0 size(tmp,2) max(find(s2m))];
            lBOX = BOX;
            lBOX(4) = lBOX(4) - 50;

            sMASTERI{e} = tmp;
            sMASTERBOX{e} = BOX;
            sMASTERBOXl{e} = lBOX;
            sqBOX{e} = {};

            ltmp = imcrop(tmp,lBOX);
            ltmp = stdfilt(ltmp,ones(11));
            ltmp = imfilter(ltmp,fspecial('average',5),'replicate');
            s1 = std(ltmp,1,1);
            s1 = bindVec(log(s1+.1));

            R = regionprops(s1 > .4,'Area','PixelIdxList');



            s1 = zeros(size(s1));
            [~,sidx] = sort([R.Area]);
            s1(R(sidx(end)).PixelIdxList) = 1;
            wididx = find(s1);
            ibk(e) = backGroundKernel(tmp);
            mBOX = [min(wididx) 0 sum(s1) lBOX(4)];
            sMASTERBOXcrop{e} = mBOX;
            jtmp = imcrop(ltmp,mBOX);
            s2 = std(jtmp,1,2);
            upper = find(s2 > 1);
            upper = min(upper);
            MBOX = mBOX;
            MBOX(2) = upper;
            MBOX(4) = BOX(4) - upper;
            sqBOX{e} = MBOX;
            skid(e) = A.kernel_id(e);
            sPOPU{e} = A.POP{e};
            %{
            if disp
                imshow(tmp,[])
                title(imageList{(e)});
                hold on

                rectangle('Position',sMASTERBOX{e},'EdgeColor','r')
                rectangle('Position',sMASTERBOXl{e},'EdgeColor','g')
                rectangle('Position',sMASTERBOXcrop{e},'EdgeColor','y')
                rectangle('Position',sqBOX{e},'EdgeColor','c')

                hold off
                drawnow
                pause(.3)
            end
            %}
        end
        
    
    
     
    e
    numel(imageList)
        
        
        
        
        
    catch ME2
        %sMASTERI{e} = [];
        %sMASTERBOX{e} = [];
    
    end
    
end
%% zip sides
close all
disp = 0;

sideDepth = [];
sideLength = [];
eidx = [];
pop_side = {};


for e = 1:numel(sMASTERI)
    if ~isempty(sMASTERI{e}) && ~isempty(sqBOX{e})
        if disp
            imshow(sMASTERI{e},[]);

            hold on

            rectangle('Position',sMASTERBOX{e},'EdgeColor','r')
            rectangle('Position',sMASTERBOXl{e},'EdgeColor','g')
            rectangle('Position',sMASTERBOXcrop{e},'EdgeColor','y')
            rectangle('Position',sqBOX{e},'EdgeColor','c')

            hold off
            drawnow
            pause(.1)
        end        
        sideDepth = [sideDepth;[sqBOX{e}(3) A.kernel_id{e}]];
        sideLength = [sideLength;[sqBOX{e}(4) A.kernel_id{e}]];
        eidx = [eidx e];
        pop_side{end+1} = A.POP{e};
    end
end
%% view sides
close all
disp = 1;


eidx = [];
kidx = [54175];
for u = 1:numel(kidx)
    eidx(u) = find(cell2mat(A.kernel_id) == kidx(u));    
end


for e = 1:numel(sMASTERI)
    if ~isempty(sMASTERI{e}) && ~isempty(sqBOX{e})
        if disp
            imshow(sMASTERI{e},[]);

            hold on

            rectangle('Position',sMASTERBOX{e},'EdgeColor','r')
            rectangle('Position',sMASTERBOXl{e},'EdgeColor','g')
            rectangle('Position',sMASTERBOXcrop{e},'EdgeColor','y')
            rectangle('Position',sqBOX{e},'EdgeColor','c')

            hold off
            drawnow
            pause(.1)
        end        
    end
end
%% zip the front
close all
disp = 0;

frontWidth = [];
frontLength = [];
pop_front = {};
IFidx = [];



for e = 1:numel(MASTERI)
    if ~isempty(MASTERI{e}) && ~isempty(qBOX{e})
        if disp
            imshow(MASTERI{e},[]);

            hold on

            rectangle('Position',MASTERBOX{e},'EdgeColor','r')
            rectangle('Position',MASTERBOXl{e},'EdgeColor','g')
            rectangle('Position',MASTERBOXcrop{e},'EdgeColor','y')
            rectangle('Position',qBOX{e},'EdgeColor','c')

            hold off
            drawnow
            pause(.1)
        end                

        frontWidth = [frontWidth;[qBOX{e}(3) A.kernel_id{e} ibk(e)]];
        frontLength = [frontLength;[qBOX{e}(4) A.kernel_id{e} ibk(e)]];
        pop_front{end+1} = A.POP{e};
        IFidx = [IFidx e];        
    end
end
%% view the front
figure;
%close all
disp = 1;

eidx = [];
kidx = [54298];
for u = 1:numel(kidx)
    eidx(u) = find(cell2mat(A.kernel_id) == kidx(u));    
end


for e = 1:numel(MASTERI)
    if ~isempty(MASTERI{e}) && ~isempty(qBOX{e})
        if disp
            imshow(MASTERI{e},[]);

            hold on

            rectangle('Position',MASTERBOX{e},'EdgeColor','r')
            rectangle('Position',MASTERBOXl{e},'EdgeColor','g')
            rectangle('Position',MASTERBOXcrop{e},'EdgeColor','y')
            rectangle('Position',qBOX{e},'EdgeColor','c')

            hold off
            drawnow
            pause(.1)
        end                   
    end
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get populate for convert to inches for NAM_parents for the front view
oldRes = strcmp(pop_front,'NAM_parents') | strcmp(pop_front,'NC-350_RILs');
oldRes = strcmp(pop_front,'NAM_parents');
kid = frontWidth(oldRes,2);
for e = 1:numel(kid)
    [num2str(frontWidth(e,1)) '-->' num2str(kid(e)) '-->' rkernelID.get(kid(e))]
end
% note kernelID's are wrong
%613-->54257-->07S-2073-05-B*_A_6--NAM_parents
%610-->54255-->07S-2073-05-B*_A_4--NAM_parents
%592-->54254-->07S-2073-05-B*_A_3--NAM_parents
%591-->54265-->07S-2073-05-B*_B_5--NAM_parents
%606-->54262-->07S-2073-05-B*_B_3--NAM_parents
%565-->54261-->07S-2073-05-B*_B_2--NAM_parents
%591-->54272-->07S-2073-05-B*_C_4--NAM_parents
%524-->54269-->07S-2073-05-B*_C_1--NAM_parents
%506-->54281-->07S-2073-05-B*_D_4--NAM_parents
%663-->54298-->07S-2073-05-B*_F_3--NAM_parents
%% make hand measurements NAM
firstImageFile = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/NAM parental lines/07S_2073_05B_07S_2072_04B_07S_2047_01B.tif';
scan1 = imread(firstImageFile);
scan1 = imcrop(scan1);
%% show NAM
figure;
imshow(scan1,[])
%% gather clicks for front view NAM
[c1 c2 V] = impixel(scan1);
WID = [c1 c2];
dW = diff(WID,1,1);
dW = sum(dW.*dW,2).^.5;
dW = dW(1:2:end);
%% make measure by hand and store in hard coded table
TABLE1 = [613;610;592;591;606;565;591;524;506;663];
TABLE1 = [TABLE1 dW];
TABLE1 = [TABLE1;[0 0]];
figure;
plot(TABLE1(:,2),TABLE1(:,1),'.')
hold on
p = polyfit(TABLE1(:,2),TABLE1(:,1),1);
pv = polyval(p,[0 1000]);
plot([0 1000],pv,'g')'
RATIO = TABLE1(:,1).*TABLE1(:,2).^-1;
p(1) = nanmean(RATIO)
p(2) = 0;
pv = polyval(p,[0 1000]);
plot([0 1000],pv,'c');
%% get populate for convert to inches for NAM_parents for the side view
oldRes = find(strcmp(pop_side,'NAM_parents'));
kid = sideDepth(oldRes,2);

for e = 1:numel(kid)
    [num2str(sideDepth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' rkernelID.get((kid(e)))]
end
% kernel ID need not be right - got them from old version of database
%519-->54257-->07S-2073-05-B*_A_6--NAM_parents
%587-->54255-->07S-2073-05-B*_A_4--NAM_parents
%609-->54254-->07S-2073-05-B*_A_3--NAM_parents
%584-->54265-->07S-2072-04-B*_B_5--NAM_parents
%605-->54262-->07S-2073-05-B*_B_3--NAM_parents
%557-->54261-->07S-2073-05-B*_B_2--NAM_parents
%472-->54272-->07S-2073-05-B*_C_4--NAM_parents
%538-->54269-->07S-2073-05-B*_C_1--NAM_parents
%% for table 2
[c12 c22 V] = impixel(scan1);
WID2 = [c12 c22];
dW2 = diff(WID2,1,1);
dW2 = sum(dW2.*dW2,2).^.5;
dW2 = dW2(1:2:end);
%% make measure by hand and store in hard coded table
TABLE2 = [519;587;609;584;605;557;472;538];
TABLE2 = [TABLE2 dW2];
%TABLE2 = [TABLE2;[0 0]];
figure;
plot(TABLE2(:,2),TABLE2(:,1),'.')
hold on
p2 = polyfit(TABLE2(:,2),TABLE2(:,1),1);
pv2 = polyval(p2,[0 1000]);
plot([0 1000],pv2,'g')'
RATIO2 = TABLE2(:,1).*TABLE2(:,2).^-1;
p2(1) = nanmean(RATIO2)
p2(2) = 0;
pv2 = polyval(p2,[0 1000]);
plot([0 1000],pv2,'c');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WIDIV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/WisconsinDiversityPanel/03May13/WISN11_22645_WISN10_13465_WISN10_13464.tif';
I = imread(fileName);
sI = imcrop(I);

%% get populate for convert to inches for back_room and therefore WIDIV for the side view
oldRes = find(strcmp(pop_side,'seedling_phenotyping_widiv'));
kid = sideDepth(oldRes,2);
tmpM = [];
for e = 1:numel(kid)
    plt = [num2str(sideDepth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' rkernelID.get((kid(e)))];
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_A_5--seedling_phenotyping_widiv')
        tmpM(1) = (sideDepth(oldRes(e),1));       
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_A_6--seedling_phenotyping_widiv')
        tmpM(2) = (sideDepth(oldRes(e),1));
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_B_5--seedling_phenotyping_widiv')
        tmpM(3) = (sideDepth(oldRes(e),1));
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_B_6--seedling_phenotyping_widiv')
        tmpM(4) = (sideDepth(oldRes(e),1));
    end
end
%% get populate for convert to inches for back_room and therefore WIDIV for the front view
oldRes = find(strcmp(pop_front,'seedling_phenotyping_widiv'));
kid = frontWidth(oldRes,2);
tmpM2 = [];
for e = 1:numel(kid)
    plt = [num2str(frontWidth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' rkernelID.get((kid(e)))];
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_A_5--seedling_phenotyping_widiv')
        tmpM2(1) = (frontWidth(oldRes(e),1));       
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_A_6--seedling_phenotyping_widiv')
        tmpM2(2) = (frontWidth(oldRes(e),1));
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_B_5--seedling_phenotyping_widiv')
        tmpM2(3) = (frontWidth(oldRes(e),1));
    end
    if strcmp(rkernelID.get((kid(e))),'WISN11_22645*_B_6--seedling_phenotyping_widiv')
        tmpM2(4) = (frontWidth(oldRes(e),1));
    end
end
%494-->79598-->WISN10_13352*_B_1--seedling_phenotyping_widiv
%538-->79587-->WISN10_13352*_B_4--seedling_phenotyping_widiv
%520-->79589-->WISN10_13352*_B_6--seedling_phenotyping_widiv
%558-->79590-->WISN10_13352*_B_7--seedling_phenotyping_widiv
%% make hand measurements NAM
secondImageFile = '/mnt/spaldingdata/Bessie/kernel scans from flatbed/WisconsinDiversityPanel/03May13/WISN10_13757_WISN10_13352_WISN11_22647.tif';
scan2 = imread(secondImageFile);
scan2 = imcrop(scan2);
%% show NAM
figure;
imshow(scan2,[])
%% gather clicks for side view
[c13 c23 V] = impixel(sI);
WID3 = [c13 c23];
dW3 = diff(WID3,1,1);
dW3 = sum(dW3.*dW3,2).^.5;
dW3 = dW3(1:2:end);
%% make measure by hand and store in hard coded table
TABLE3 = tmpM';
TABLE3 = [TABLE3 dW3];
figure;
plot(TABLE3(:,2),TABLE3(:,1),'.')
hold on
p3 = polyfit(TABLE3(:,2),TABLE3(:,1),1);
pv3 = polyval(p3,[0 1000]);
plot([0 1000],pv3,'g')'
RATIO3 = TABLE3(:,1).*TABLE3(:,2).^-1;
p3(1) = nanmean(RATIO3)
p3(2) = 0;
pv3 = polyval(p3,[0 1000]);
plot([0 1000],pv3,'c');
%% gather clicks for front view
[c14 c24 V] = impixel(sI);
WID4 = [c14 c24];
dW4 = diff(WID4,1,1);
dW4 = sum(dW4.*dW4,2).^.5;
dW4 = dW4(1:2:end);
%% make measure by hand and store in hard coded table
close all
TABLE4 = tmpM2';
TABLE4 = [TABLE4 dW4];
figure;
plot(TABLE4(:,2),TABLE4(:,1),'.')
hold on
p4 = polyfit(TABLE4(:,2),TABLE4(:,1),1);
pv4 = polyval(p4,[0 1000]);
plot([0 1000],pv4,'g')'
RATIO4 = TABLE4(:,1).*TABLE4(:,2).^-1;
p4(1) = nanmean(RATIO4)
p4(2) = 0;
pv4 = polyval(p4,[0 1000]);
plot([0 1000],pv4,'c');

