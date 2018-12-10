%% notes from 2014.09.26
% refocus on NAM lines
% refocus on shape data
%% try to position to well names
k = translateWellNames(f);
%% load NAM from database
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...
                'JOIN predictions ' ...
                'on kernels.id = predictions.kernel_id ' ...                
                'JOIN averageweightspectra_vw ' ...
                'ON kernels.id=averageweightspectra_vw.kernel_id ' ...
                'JOIN kernel_3d ' ...
                'ON kernels.id=kernel_3d.kernel_id ' ...                
                'JOIN reporting.kernel_dims_report_tbl '...
                'on kernels.id=reporting.kernel_dims_report_tbl.kernel_id '...
                'JOIN files ' ...
                'on kernels.img_gray_side_fid = files.id ' ...                
                'WHERE population_lines.type = ' '''NAM_parents'' '];
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
%% split data NAM
fidx = [];
fidx(1) = find(strcmp(fieldString,'top_mi_area'));
fidx(2) = find(strcmp(fieldString,'top_major_axis'));
fidx(3) = find(strcmp(fieldString,'top_minor_axis'));

fidx(4) = find(strcmp(fieldString,'front_mi_area'));
fidx(5) = find(strcmp(fieldString,'front_major_axis'));
fidx(6) = find(strcmp(fieldString,'front_minor_axis'));

fidx(7) = find(strcmp(fieldString,'side_mi_area'));
fidx(8) = find(strcmp(fieldString,'side_major_axis'));
fidx(9) = find(strcmp(fieldString,'side_minor_axis'));


%f.specData = cell2mat(results(:,54:54+781));
f.specData = cell2mat(results(:,(54-3):54+781));
f.genoType = results(:,16);
f.kernel_id = cell2mat(results(:,34));
f.predictions = cell2mat(results(:,[35:40 836]));
f.specData = f.specData(:,4:end);
f.shapeData = cell2mat(results(:,fidx));
f.plateName = results(:,6);
f.position = results(:,23);
f.shapeDataL = cell2mat(results(:,881:884));
f.xpos = cell2mat(results(:,24));
f.imageFRONT = results(:,end);
f.imageSIDE = results(:,end);
f.imageTOP = results(:,end);
f.pos_x = cell2mat(results(:,find(strcmp(fieldString,'cob_position_x'))));
%% clean data NAM
GT = 'Hp301';
%GT = '';
ridx = find(any(isnan(f.specData),2) | ...            
            any(f.shapeData(:,4) > 3*10^5,2) | ...
            any(f.shapeData(:,5) > 400,2)  | ...
            any(f.shapeDataL(:,end) > .5*10^9) | ...
            f.predictions(:,5) > f.predictions(:,6) | ... 
            isnan(f.xpos) | ... 
            any(strcmp('W22^ACR',f.genoType),2) | ...
            any(strcmp(GT,f.genoType),2) ...
            ); 
f.specData(ridx,:) = [];
f.kernel_id(ridx) = [];
f.genoType(ridx) = [];
f.predictions(ridx,:) = [];
f.shapeData(ridx,:) = [];
f.plateName(ridx,:) = [];
f.position(ridx,:) = [];
f.shapeDataL(ridx,:) = [];
f.xpos(ridx,:) = [];
f.imageTOP(ridx) = [];
f.pos_x(ridx) = [];
%% labels
%LAB.shapeOLD = fieldString(1005:1008);
LAB.shape = {'Height' 'Width' 'Depth' 'Volume'};
LAB.spec = [fieldString([35:40 836]);'mgOil'];
LAB.spec = {'Percent Protein', 'Percent Oil','Total Density','Material Density','Material Volume','Total Volume','Weight','mg Oil'};
LAB.gravi = {'Swing' 'Final'};
LAB.length = {'Length' 'Final'};
LAB.shapePCA = {'seed set','kernel width','seed size'};
%% add air space
airSPACE = f.predictions(:,6) - f.predictions(:,5);
f.predictions = [f.predictions airSPACE];
LAB.spec{end+1} = 'airspace';
%% remove air space
f.predictions(:,end) = [];
LAB.spec(end) = [];
%% add area and ratio to shape
f.shapeDataL = [f.shapeDataL f.shapeDataL(:,2).*f.shapeDataL(:,3) f.shapeDataL(:,2).*f.shapeDataL(:,3).^-1];
LAB.shape{end+1} = 'cap area';
LAB.shape{end+1} = 'cap ratio';
%% remove area and ratio to shape
f.shapeDataL(:,end-1:end) = [];
LAB.shape{end} = [];
LAB.shape{end} = [];
%% analysis for FIGURE 2 for paper
T = f;
close all
X = [T.specData];
Y = [T.shapeDataL(:,[1:3])];

Y = zscore(Y);

[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
%myCCA(X,Y)
% coeffs for X
for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end
% coeffs for Y
for e = 1:size(mU,2)
    figure;
    bar(mB(:,e));
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end

% struct core for Y
[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.shape(1:3),'XTick',1:numel(LAB.shape(1:3)))
    title(['Core Shape:' num2str(e)]);
end
% struct core for X
[Xcore XcoreP] = corr(mU,X);

for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end
%% create correlation figures for 3 biplot 1-2
xDIM = 1;
yDIM = 2;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';


figure;
hold on
C = [c1' c2'];


LEG = {};
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,1:2);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        %[tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        %plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
        LEG{end+1} = UQ{u};
    end
end


plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');


% try to double check
for e = 1:size(mU,2)
    for l = 1:size(f.predictions,2)
        CHECK(e,l)= corr(f.predictions(:,l),mU(:,e));
    end
end
CHECK = corr(f.predictions,C);
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,1)
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end


% try to double check
CHECK = [];

for e = 1:size(mU,2)
    for l = 1:size(Y,2)
        CHECK(e,l)= corr(Y(:,l),mV(:,e));
    end
end
CHECK = corr(Y,C);
L = {'Height' 'Width' 'Depth' 'Cap Area' 'Cap Ratio'};
C_lessP = CHECK;
for e = 1:size(C_lessP,1)    
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
legend(LEG)

axis([-6 6 -6 6])
%% create correlation figures for 3 biplot 2-3

xDIM = 2;
yDIM = 3;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';


figure;
hold on
C = [c1' c2'];

% raster the genotypes
LEG = {};
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,1:2);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        %[tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        %plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
        LEG{end+1} = UQ{u};
    end
end
%{
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,1:2);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        [tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        plot(tx,ty,'k--','LineWidth',.2);                
    end
end
%}


plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');

% biplot correlations
CHECK = corr(f.predictions,C);
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,1)
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end

% biplot correlations
CHECK = corr(Y,C);
L = {'Height' 'Width' 'Depth'};
C_lessP = CHECK;
for e = 1:size(C_lessP,1)    
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
legend(LEG)

axis([-6 6 -6 6])
%% create correlation figures for 3 biplot 2-3

xDIM = 1;
yDIM = 2;
zDIM = 3;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
c3 = [2^-.5 2^-.5]*[mU(:,zDIM) mV(:,zDIM)]';


figure;
hold on
C = [c1' c2' c3'];

% raster the genotypes
LEG = {};
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,:);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        %[tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        %plot(tx,ty,'k--','LineWidth',.2);
        plot3(mean(sub(:,1)),mean(sub(:,2)),mean(sub(:,3)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
        LEG{end+1} = UQ{u};
    end
end


plot3(C(:,1),C(:,2),C(:,3),'.','MarkerSize',.2,'Color','k');

% biplot correlations
CHECK = corr(f.predictions,C);
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,1)
    P = C_lessP(e,:);
    quiver3(0,0,0,P(1),P(2),P(3),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    tZ = linspace(0,textMAG*P(3),100);
    plot3(tX,tY,tZ,'--','LineWidth',1,'Color','k');
end

% biplot correlations
CHECK = corr(Y,C);
L = {'Height' 'Width' 'Depth'};
C_lessP = CHECK;
for e = 1:size(C_lessP,1)    
    P = C_lessP(e,:);
    quiver3(0,0,0,P(1),P(2),P(3),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    tZ = linspace(0,textMAG*P(3),100);
    plot3(tX,tY,tZ,'--','LineWidth',1,'Color','r');
end
legend(LEG)

axis([-6 6 -6 6 -6 6])
%% load nearest n
NN = readtext('/mnt/spaldingdata/nate/NearestNeighbors.csv');
NN(1,:) = [];
n = [];
for e = 1:size(NN,1)
    if ~isempty(str2num(NN{e,6}(2:end-1)))
        n = [n;[str2num(NN{e,5}(2:end-1)) str2num(NN{e,6}(2:end-1))]];
    end
end
%% cap area with cvN
xDIM = 1;
yDIM = 2;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
C = [c1' c2'];

CAPArea = f.shapeDataL(:,2).*f.shapeDataL(:,3);
corr(CAPArea,C(:,2))

xDIM = 2;
yDIM = 3;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
C = [c1' c2'];
CAPArea = f.shapeDataL(:,2).*f.shapeDataL(:,3).^-1;
corr(CAPArea,C(:,2))
%% plot over ears
xDIM = 1;
yDIM = 2;
zDIM = 3; 
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
c3 = [2^-.5 2^-.5]*[mU(:,zDIM) mV(:,zDIM)]';
C = [c1' c2' c3'];
tmp = {};
UQ = unique(f.genoType);
for e = 1:numel(UQ)    
    fidx = find(strcmp(f.genoType,UQ{e}));    
    [XPOS sidx] = sort(f.pos_x(fidx));
    tmp{e} = C(fidx,:);
    tmp{e} = f.predictions(fidx,4);
    tmp{e} = tmp{e}(sidx,:);
    [XPOS selidx] = unique(XPOS);
    %YC(:,:,e) = interp1(XPOS,tmp{e}(selidx,:),linspace(min(XPOS),max(XPOS),100),'spline');
    plot(XPOS,tmp{e}(selidx,:))
    title(UQ{e});
    waitforbuttonpress
    tmpy = tmp{e};    
    SD(e,:) = std(tmpy,1,1);
    U(e,:) = mean(tmpy,1);
    nSD(e,:) = SD(e,:).*U(e,:).^-1;
    pause(.1);
    drawnow
end
%% plot max var
%% plot min var
[J midx] = max(SD,[],1);
for e = 1:numel(midx)
    figure
    plot(tmp{midx(e)})
    title([UQ{midx(e)} 'max-' num2str(e)])
    legend({'cv1','cv2','cv3'});
end
%% plot min var
[J midx] = min(SD,[],1);
for e = 1:numel(midx)
    figure
    plot(tmp{midx(e)})
    title([UQ{midx(e)} 'min-' num2str(e)])
    legend({'cv1','cv2','cv3'});
end
%% plot max mean
[J midx] = max(U,[],1);
for e = 1:numel(midx)
    figure
    plot(tmp{midx(e)})
    title([UQ{midx(e)} 'max-' num2str(e)])
    legend({'cv1','cv2','cv3'});
end
%% plot min  mean
[J midx] = min(U,[],1);
for e = 1:numel(midx)
    figure
    plot(tmp{midx(e)})
    title([UQ{midx(e)} 'min-' num2str(e)])
    legend({'cv1','cv2','cv3'});
end


%% analysis of nearest and stacking
newMeasure = [];
for e = 1:size(n,1)
    fidx = find(f.kernel_id == n(e,1));
    if ~isempty(fidx)
        %newMeasure = [newMeasure ;[C(fidx,:) f.predictions(fidx,:) f.shapeDataL(fidx,:) n(e,2)]];
        vol = f.predictions(fidx,6);
        newMeasure = [newMeasure ;[C(fidx,:) f.shapeDataL(fidx,2)*f.shapeDataL(fidx,3) f.shapeDataL(fidx,2)/f.shapeDataL(fidx,3) f.shapeDataL(fidx,1) n(e,2)]];
    end
end
corr(newMeasure)
%% view from 3 axis hand measurement
IDX = 1987;
for IDX = 1:10
            try
                ext = '_F.tiff';            
                tmp = imread([f.imageTOP{IDX} ext]);                
                ext = '_S2.tiff';
                tmp2 = imread([f.imageTOP{IDX} ext]);
                ext = '_T.tiff';
                tmp3 = imread([f.imageTOP{IDX} ext]);
            catch ME2
                try
                    ext = '_F.tiff';            
                    tmp = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));
                    ext = '_S2.tiff';
                    tmp2 = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));
                    ext = '_T.tiff';
                    tmp3 = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));
                catch ME2
                    try
                        ext = '_F.tiff';            
                        tmp = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);
                        ext = '_S2.tiff';
                        tmp2 = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);
                        ext = '_T.tiff';
                        tmp3 = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);
                    catch ME2
                        ext = '_F.tiff';            
                        tmp = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));
                        ext = '_S2.tiff';
                        tmp2 = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));
                        ext = '_T.tiff';
                        tmp3 = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));
                        
                    end
                end
            end

figure;imshow(tmp,[])
[p1 p2 V] = impixel();
D1 = p1(2) - p2(2);
figure;imshow(tmp2,[])
figure;imshow(tmp3,[])
[p1 p2 V] = impixel();
D2 = p1(1) - p2(1);
DT(IDX) = abs(D2 - D1);
end
%% rank kernels and view images
UQ = unique(f.genoType);
rankVec = f.shapeDataL(:,3);
rankVec = C(:,2);
secondaryIndex = f.shapeDataL;
uVec = [];uVec2 = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    uVec(u) = mean(rankVec(fidx));
    uVec2(u,:) = mean(secondaryIndex(fidx,:),1);
end
[J sidx] = sort(uVec,'descend');
uVec = uVec(sidx);
uVec2 = uVec2(sidx,:);
UQ = UQ(sidx);
tI = [];
lI = [];BOX = [];ext = '_S2.tiff';
d = figure;
for u = [1:4 numel(UQ)-3:numel(UQ)]    
%for u = [1:numel(UQ)]    
    fidx = find(strcmp(UQ{u},f.genoType));
    distance = abs(rankVec(fidx) - uVec(u));
    [distance sidx] = sort(distance);
    
    %r = randperm(numel(fidx));
    %fidx = fidx(r);
    
    sI = [];cnt = 1;e = 1;
    while cnt ~= 2
        IDX = fidx(sidx(1));
        try
            try
                tmp = imread([f.imageTOP{IDX} ext]);                               
            catch ME2
                try
                    tmp = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));              
                catch ME2
                    try
                        tmp = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);                
                    catch ME2
                        tmp = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));                        
                    end
                end
            end
            
            if strcmp(ext,'_S2.tiff')
                tmp(end-245:end,:) = [];    
            else
                tmp(end-100:end,:) = [];    
            end
            
            
            
            E = edge(tmp,'canny',[],5);
            E = bwareaopen(E,10);
            [H,T,R] = hough(E,'Theta',-90);
            P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:))))            
            lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
            imshow(E,[]);
            hold on
            max_len = 0;
            for k = 1:length(lines)
               xy = [lines(k).point1; lines(k).point2];               
               
               plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

               % Plot beginnings and ends of lines
               plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
               plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            
               
               % Determine the endpoints of the longest line segment
               len = norm(lines(k).point1 - lines(k).point2);
               if ( len > max_len)
                  max_len = len;
                  xy_long = xy;
               end
            end
            BOX = [0 0 size(tmp,2) xy_long(3)];
            
            
            
            level = graythresh(tmp);
            btmp = double(tmp)/255 > level;
            
            if isempty(BOX)
                figure;
                [J BOX] = imcrop(tmp);
                close all;
            else
                figure(d);
            end
            hold off
            tmp = imcrop(tmp,BOX);            
            imshow(tmp,[]);
            
            PER = .2;
            u1 = mean(tmp,1);            
            s1 = std(double(tmp),1,1);
            s1 = bindVec(s1);
            level = graythresh(s1);
            level = level - PER*level;
            idx1 = find(s1 > level);
            ux = mean(idx1);
            MAX1 = max(idx1);
            MIN1 = min(idx1);
            
            
            u2 = mean(tmp,2);            
            s2 = std(double(tmp),1,2);
            s2 = bindVec(s2);
            level = graythresh(s2);
            level = level - PER*level;
            idx2 = find(s2 > level);
            ux = mean(idx2);
            MAX2 = max(idx2);
            MIN2 = min(idx2);
            MAX2 = size(tmp,1);
            
            
            RECTVALUE = f.shapeDataL(IDX,[1 3]);
            
            rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
            BOTLEFT = [MIN1 MAX2];
            UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
            rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');
            
            nBOX = [MIN1,MIN2,MAX1-MIN1,MAX2-MIN2];
            img = imcrop(tmp,nBOX);
            
            sIL = padarray(img,[SZ(1)-size(img,1) 0],'pre');
            sIT = padarray(img,[0 SZ(2)-size(img,2)],'post');
            %tmp = padarray(tmp,[0 SZ(2)-size(tmp,2)],'post');
            drawnow;            
            %sI = [sI tmp];
            e = e + 1;cnt = cnt +1;
            %waitforbuttonpress
            hold off
        catch ME;e = e + 1;
            ME;
        end
    end
    %I = [I;sI];
    lI = [lI,sIL];
    tI = [tI;sIT];
end
imshow(lI,[]);hold on
YVEC = uVec2(:,3);
plot([1 size(lI,2)],[size(lI,1)-YVEC(1) size(lI,1)-YVEC(end)],'g');
plot([1 size(lI,2)],[size(lI,1)-YVEC(1) size(lI,1)-YVEC(1)],'r');
plot([1 size(lI,2)],[size(lI,1)-YVEC(end) size(lI,1)-YVEC(end)],'b');
figure;
imshow(tI,[]);
%% rank kernels and view images -- hold constant
fidx = find(abs(C(:,1)) < .2);
[JUNK sidx] = sort(C(fidx,2),'descend');
fidx = fidx(sidx);
fidx = fidx(1:20:end);
tI = [];
lI = [];BOX = [];ext = '_S2.tiff';
d = figure;
PER = .6;
%PER = .2;
COMP = [];
for i = 1:numel(fidx);
    sI = [];cnt = 1;e = 1;
    while cnt ~= 2
        IDX = fidx(i);
        COMP = [COMP f.predictions(IDX,2)];
        try
            try
                tmp = imread([f.imageTOP{IDX} ext]);                               
            catch ME2
                try
                    tmp = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));              
                catch ME2
                    try
                        tmp = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);                
                    catch ME2
                        tmp = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));                        
                    end
                end
            end
            
            
            SZ = size(tmp);
            
            if strcmp(ext,'_S2.tiff')
                tmp(end-245:end,:) = [];    
            else
                tmp(end-100:end,:) = [];    
            end
            
            
            
            E = edge(tmp,'canny');
            E = bwareaopen(E,10);
            [H,T,R] = hough(E,'Theta',-90);
            P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:))))            
            lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
            imshow(E,[]);
            
            
            % 
            hold on
            max_len = 0;
            for k = 1:length(lines)
               xy = [lines(k).point1; lines(k).point2];               
               
               plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

               % Plot beginnings and ends of lines
               plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
               plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            
               
               % Determine the endpoints of the longest line segment
               len = norm(lines(k).point1 - lines(k).point2);
               if ( len > max_len)
                  max_len = len;
                  xy_long = xy;
               end
            end
            BOX = [0 0 size(tmp,2) xy_long(3)];
            
            
            
            level = graythresh(tmp);
            btmp = double(tmp)/255 > level;
            
            if isempty(BOX)
                figure;
                [J BOX] = imcrop(tmp);
                close all;
            else
                figure(d);
            end
            hold off
            tmp = imcrop(tmp,BOX);            
            
            imshow(tmp,[]);
            
            %PER = .9;
            u1 = mean(tmp,1);            
            s1 = std(double(tmp),1,1);
            s1 = bindVec(s1);
            level = graythresh(s1);
            level = level - PER*level;
            idx1 = find(s1 > level);
            ux = mean(idx1);
            MAX1 = max(idx1);
            MIN1 = min(idx1);
            
            
            u2 = mean(tmp,2);            
            s2 = std(double(tmp),1,2);
            s2 = bindVec(s2);
            level = graythresh(s2);
            level = level - PER*level;
            idx2 = find(s2 > level);
            ux = mean(idx2);
            MAX2 = max(idx2);
            MIN2 = min(idx2);
            MAX2 = size(tmp,1);
            
            
            RECTVALUE = f.shapeDataL(IDX,[2 3]);
            
            rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
            BOTLEFT = [MIN1 MAX2];
            UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
            rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');
            
            nBOX = [MIN1,MIN2,MAX1-MIN1,MAX2-MIN2];
            img = imcrop(tmp,nBOX);
            %hello = imrotate(img,-90);
            hello = img;
            sIL = padarray(hello,[SZ(1)-size(hello,1) 0],'pre');
            sIT = padarray(img,[0 SZ(2)-size(img,2)],'post');
            %tmp = padarray(tmp,[0 SZ(2)-size(tmp,2)],'post');
            drawnow;            
            %sI = [sI tmp];
            e = e + 1;cnt = cnt +1;
            %waitforbuttonpress
            hold off
        catch ME;e = e + 1;
            ME;
        end
    end
    %I = [I;sI];
    lI = [lI,sIL];
    tI = [tI;sIT];
    
end
figure;
imshow(lI,[]);hold on
%imwrite(lI,'/mnt/spaldingdata/nate/communications/papers/paper with jeffG/figure1/kernelAlongCv1.tif');
imwrite(lI,'/mnt/spaldingdata/nate/communications/papers/paper with jeffG/figure1/kernelAlongCv2.tif');
%imwrite(tI,'/mnt/spaldingdata/nate/communications/papers/paper with jeffG/figure1/kernelAlongCv3.tif');
%YVEC = uVec2(:,3);
%plot([1 size(lI,2)],[size(lI,1)-YVEC(1) size(lI,1)-YVEC(end)],'g');
%plot([1 size(lI,2)],[size(lI,1)-YVEC(1) size(lI,1)-YVEC(1)],'r');
%plot([1 size(lI,2)],[size(lI,1)-YVEC(end) size(lI,1)-YVEC(end)],'b');
figure;
imshow(tI,[]);
%% rank kernels along composition
fidx = find(abs(C(:,1)) < .5);
[JUNK sidx] = sort(C(fidx,2),'descend');
fidx = fidx(sidx);
fidx = fidx(1:1:end);
COMP = [];
for i = 1:numel(fidx);
    sI = [];cnt = 1;e = 1;    
    IDX = fidx(i);
    COMP = [COMP;[C(IDX,2) f.predictions(IDX,2)]];
end
%% create coll
UQ = unique(f.genoType);
rankVec = f.shapeDataL(:,3);
rankVec = C(:,2);
secondaryIndex = f.shapeDataL;
uVec = [];uVec2 = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    uVec(u) = mean(rankVec(fidx));
    uVec2(u,:) = mean(secondaryIndex(fidx,:),1);
end
[J sidx] = sort(uVec);
uVec = uVec(sidx);
uVec2 = uVec2(sidx,:);
UQ = UQ(sidx);
I = [];BOX = [];ext = '_S2.tiff';
d = figure;

%for u = [1:4 numel(UQ)-3:numel(UQ)]    
for u = [1:numel(UQ)]    
    fidx = find(strcmp(UQ{u},f.genoType));
    distance = abs(rankVec(fidx) - uVec(u));
    [distance sidx] = sort(distance);
    
    r = randperm(numel(fidx));
    fidx = fidx(r);
    
    sI = [];cnt = 1;e = 1;
    while cnt ~= 10
        IDX = fidx(e);
        try
            try
                tmp = imread([f.imageTOP{IDX} ext]);                                
            catch ME2
                try
                    tmp = imread(strrep([f.imageTOP{IDX}],'_S2.tiff',ext));                    
                catch ME2
                    try
                        tmp = imread([strrep(f.imageTOP{IDX},'gray/07S','gray/07') ext]);                        
                    catch ME2
                        tmp = imread(strrep([strrep(f.imageTOP{IDX},'gray/07S','gray/07')],'_S2.tiff',ext));                        
                    end
                end
            end
            
            if strcmp(ext,'_S2.tiff')
                tmp(end-245:end,:) = [];    
            else
                tmp(end-100:end,:) = [];    
            end
            %{
            
            
            E = edge(tmp,'canny',[],5);
            E = bwareaopen(E,10);
            [H,T,R] = hough(E,'Theta',-90);
            P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:))))            
            lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
            imshow(E,[]);
            hold on
            max_len = 0;
            for k = 1:length(lines)
               xy = [lines(k).point1; lines(k).point2];               
               
               plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

               % Plot beginnings and ends of lines
               plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
               plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            
               
               % Determine the endpoints of the longest line segment
               len = norm(lines(k).point1 - lines(k).point2);
               if ( len > max_len)
                  max_len = len;
                  xy_long = xy;
               end
            end
            BOX = [0 0 size(tmp,2) xy_long(3)];
            %}
            %waitforbuttonpress
            
            
            
            level = graythresh(tmp);
            btmp = double(tmp)/255 > level;
            
            if isempty(BOX)
                figure;
                [J BOX] = imcrop(tmp);
                close all;
            else
                figure(d);
            end
            hold off
            tmp = imcrop(tmp,BOX);            
            imshow(tmp,[]);
            
            PER = .2;
            u1 = mean(tmp,1);            
            s1 = std(double(tmp),1,1);
            s1 = bindVec(s1);
            level = graythresh(s1);
            level = level - PER*level;
            idx1 = find(s1 > level);
            ux = mean(idx1);
            MAX1 = max(idx1);
            MIN1 = min(idx1);
            
            
            u2 = mean(tmp,2);            
            s2 = std(double(tmp),1,2);
            s2 = bindVec(s2);
            level = graythresh(s2);
            level = level - PER*level;
            idx2 = find(s2 > level);
            ux = mean(idx2);
            MAX2 = max(idx2);
            MIN2 = min(idx2);
            MAX2 = size(tmp,1);
            
            
            RECTVALUE = f.shapeDataL(IDX,[1 3]);
            
            rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
            BOTLEFT = [MIN1 MAX2];
            UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
            rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');
            
            nBOX = [MIN1,MIN2,MAX1-MIN1,MAX2-MIN2];
            img = imcrop(tmp,nBOX);
            
            img = padarray(img,[SZ(1)-size(img,1) 0],'pre');
            img = padarray(img,[0 SZ(2)-size(img,2)],'post');
            %tmp = padarray(tmp,[0 SZ(2)-size(tmp,2)],'post');
            drawnow;            
            sI = [sI tmp];
            e = e + 1;cnt = cnt +1;
            %waitforbuttonpress
            hold off
        catch ME;e = e + 1;
            ME;
        end
    end
    I = [I;sI];
    %lI = [lI,sIL];
    %tI = [tI;sIT];
end
imshow(I,[])


%% analysis on all images from side
ext = '_S2.tiff';
Pnew = zeros(numel(f.imageTOP),2);
imageList = A.imageTOP;
clear tmp;

for e = 1:numel(imageList)
    try
        try
            tmp = imread([imageList{(e)} ext]);
        catch ME2
            try
                tmp = imread([imageList{(e)}]);
            catch ME2
                try
                    tmp = imread([strrep(imageList{(e)},'gray/07S','gray/07') ext]);
                catch ME2
                    tmp = imread([strrep(imageList{(e)},'gray/07S','gray/07')]);
                end
            end
        end
    catch ME2
    
    end
    if backGroundKernel(tmp) > 255/2 
        tmp = imcomplement(tmp);
    end
    
    if strcmp(ext,'_S2.tiff')
        tmp(end-245:end,:) = [];    
    else
        tmp(end-100:end,:) = [];    
    end
    E = edge(tmp,'canny');
    E = bwareaopen(E,10);
    [H,T,R] = hough(E,'Theta',-90);
    P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:)))) ;           
    lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
    %imshow(E,[]);
    %hold on
    max_len = 0;
    xy_long = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];               

       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    BOX = [0 0 size(tmp,2) xy_long(3)];
        

    
    level = graythresh(tmp);
    btmp = double(tmp)/255 > level;

    if isempty(BOX)
        figure;
        [J BOX] = imcrop(tmp);
        close all;
    else
        %figure(d);
    end

    tmp = imcrop(tmp,BOX);
    %imshow(tmp,[])
    
    PER = .2;
    u1 = mean(tmp,1);            
    s1 = std(double(tmp),1,1);
    s1 = bindVec(s1);
    level = graythresh(s1);
    level = level - PER*level;
    idx1 = find(s1 > level);
    ux = mean(idx1);
    MAX1 = max(idx1);
    MIN1 = min(idx1);


    u2 = mean(tmp,2);            
    s2 = std(double(tmp),1,2);
    s2 = bindVec(s2);
    level = graythresh(s2);
    level = level - PER*level;
    idx2 = find(s2 > level);
    ux = mean(idx2);
    MAX2 = max(idx2);
    MIN2 = min(idx2);
    MAX2 = size(tmp,1);


    %RECTVALUE = f.shapeDataL((e),[1 3]);

    %rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
    %BOTLEFT = [MIN1 MAX2];
    %UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
    %rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');



    %Pnew = [Pnew;[MAX1-MIN1,MAX2-MIN2]];
    %Pold = [Pold;f.shapeDataL((e),[1 3])];
    Pnew(e,:) = [MAX1-MIN1,MAX2-MIN2];    
    %Pold(e,:) = [f.shapeDataL((e),[1 3])];    
    %drawnow;            
    %sI = [sI tmp];
    %e = e + 1;cnt = cnt +1;
    e
    numel(imageList)
end
%% analysis on all images from front
close all
ext = '_F.tiff';
PnewF = zeros(numel(f.imageTOP),2);
imageList = f.imageTOP;
imageList = imageList(randperm(numel(imageList)));
clear tmp;
PER = .65;
disp = 1;

UQ = unique(f.genoType);
widx = [];
for u = 1:numel(UQ)
    cnt = 1;
    for e = 1:numel(f.genoType)
        
        if strcmp(f.genoType{e},UQ{u}) & cnt < 5
            widx = [widx e];
            cnt = cnt + 1;
        end        
    end
end
%parfor e = 1:numel(f.imageTOP)
for e = 1:numel(imageList)
    try
        try
            tmp = imread([imageList{(e)} ext]);
        catch ME2
            try
                tmp = imread(strrep([imageList{(e)}],'_S2.tiff',ext));
            catch ME2
                try
                    tmp = imread([strrep(imageList{(e)},'gray/07S','gray/07') ext]);
                catch ME2
                    tmp = imread(strrep([strrep(imageList{(e)},'gray/07S','gray/07')],'_S2.tiff',ext));
                end
            end
        end
    catch ME2
    
    end
    
    
    if backGroundKernel(tmp) > 255/2 
        tmp = imcomplement(tmp);
    end
    
    
    tmp(end-100:end,:) = [];
    E = edge(tmp,'canny');
    E = bwareaopen(E,10);
    [H,T,R] = hough(E,'Theta',-90);
    P = houghpeaks(H,1,'threshold',ceil(0.3*max(H(:)))) ;           
    lines = houghlines(E,T,R,P,'FillGap',200,'MinLength',50);
    %imshow(E,[]);
    %hold on
    max_len = 0;
    xy_long = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];               

       %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       %plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       %plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');            

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    BOX = [0 0 size(tmp,2) xy_long(3)];
        

    level = graythresh(tmp);
    btmp = double(tmp)/255 > level;

    if isempty(BOX)
        figure;
        [J BOX] = imcrop(tmp);
        close all;
    else
        %figure(d);
    end

    
    % modif crop BOX by hand
    BOX(4) = BOX(4) - 50;
    tmp = imcrop(tmp,BOX);
    
    
    tmp1 = stdfilt(tmp,ones(11));
    %tmp1 = imfilter(tmp1,fspecial('average',5));
    u1 = mean(tmp1,1);            
    s1 = std(double(tmp1),1,1);
    u1 = bindVec(u1);
    s1 = bindVec(s1);
    %s1 = s1.*u1;    
    s1 = u1;
    s1 = bindVec(s1)+.1;
    s1 = log(s1);
    s1 = bindVec(s1);
    level = graythresh(s1);
    level = level - PER*level;
    level = .10;
    binaryL = s1 > level;
    R = regionprops(binaryL,'Area','PixelIdxList');
    [~,sidx] = sort([R.Area]);
    myThresh = zeros(size(binaryL));
    myThresh(R(sidx(end)).PixelIdxList) = 1;
    
    idx1 = find(myThresh);
    
    ux = mean(idx1);
    MAX1 = max(idx1);
    MIN1 = min(idx1);

    u2 = mean(tmp,2);            
    s2 = std(double(tmp),1,2);
    s2 = bindVec(s2);
    level = graythresh(s2);
    level = level - PER*level;
    idx2 = find(s2 > level);
    ux = mean(idx2);
    MAX2 = max(idx2);
    MIN2 = min(idx2);
    MAX2 = size(tmp,1);


    %RECTVALUE = f.shapeDataL((e),[1 3]);
    if disp
        imshow(tmp1)    
        rectangle('Position',[MIN1,MIN2,MAX1-MIN1,MAX2-MIN2],'EdgeColor','r');
        drawnow
    end
    %BOTLEFT = [MIN1 MAX2];
    %UPPERLEFT = [MIN1 MAX2 - RECTVALUE(2)];
    %rectangle('Position',[UPPERLEFT(1),UPPERLEFT(2),RECTVALUE(1),RECTVALUE(2)],'EdgeColor','g');



    %Pnew = [Pnew;[MAX1-MIN1,MAX2-MIN2]];
    %Pold = [Pold;f.shapeDataL((e),[1 3])];
    PnewF(e,:) = [MAX1-MIN1,MAX2-MIN2];   
    %Pold(e,:) = [f.shapeDataL((e),[1 3])];    
    %drawnow;            
    %sI = [sI tmp];
    %e = e + 1;cnt = cnt +1;
    e
    numel(imageList)
end
PnewF(:,2) = PnewF(:,2) + 50;
%%
p = polyfit(Pnew(:,2),PnewF(:,2),1);
y = polyval(p,linspace(min(Pnew(:,2)),max(Pnew(:,2)),100));
%% over assign the shape vectors
f.shapeDataL = [Pnew(:,1) PnewF(:,1) mean([Pnew(:,2) PnewF(:,2)],2)];
f.shapeDataL = [Pnew(:,1) PnewF(:,1) mean([Pnew(:,2) Pnew(:,2)],2)];
%% check for simple examples
fidx = find(strcmp('B73',f.genoType));
figure;
plot(f.shapeDataL(:,1),f.shapeDataL(:,2),'b.');
hold on;
plot(f.shapeDataL(fidx,1),f.shapeDataL(fidx,2),'r.');

%% perform calibration for oil, etc
[D] = readtext('/mnt/spaldingdata/nate/121213_NIRCalibration_Spectra_NIRTube2_5reps_edited_AVG.csv');
X = cell2mat(D(2:end,26:end));
Y = D(2:end,[15:25]);
toRM = [];
for e = 1:size(Y,1)
    toRM(e) = 0;
    for j = 1:size(Y,2)
        toRM(e) = toRM(e) | isempty(Y{e,j});
    end
end
ridx = find(toRM);
Y(ridx,:) = [];
X(ridx,:) = [];
Y = cell2mat(Y);
%Y(:,7) = Y(:,1).*Y(:,9).^-1;
BETA = [];

for e = 1:size(Y,2)
    [XL,YL,XS,YS,BETA{e},PCTVAR,MSE] = plsregress(X,Y(:,e),8);
end
%predictVectors = PCA_BKPROJ(BETA(2:end,:)',caliE,caliU);
predictVectors = BETA(1:end,:)';
BETA = cell2mat(BETA);
preTest = [ones(size(f.specData,1),1) f.specData(:,4:end)]*BETA;
preTest2 = [ones(size(A.specData,1),1) A.specData(:,4:end)]*BETA;
NMS = D(1,[15:25]);
%preTest = [ones(size(X,1),1) X]*BETA;
%% build and add starch
sd = readtext('/home/nate/Downloads/Starch-NIRdata.csv');
%% add starch
sY = cell2mat(sd(2:361,3:4));
sS = cell2mat(sd(2:361,6:end));
sS(:,end) = [];
[XL,YL,XS,YS,beta,PCTVAR] = plsregress(sS,sY,15);
%sP = [ones(size(sS,1),1) sS]*beta;
preTest = [preTest [ones(size(f.specData,1),1) f.specData(:,4:end)]*beta];
NMS{end+1} = 'perS';
NMS{end+1} = 'mgS';
%% shotgunA
close all
figure;
%SHOT = [mBCC(:,9).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,10).*(mBCC(:,5)-mBCC(:,4)),mBCC(:,11).*(mBCC(:,5)-mBCC(:,4))];
%SHOT = [mBCC(:,1),mBCC(:,2),mBCC(:,3)];
SHOT = [preTest(:,5)/100.*preTest(:,1) preTest(:,4)/100.*preTest(:,1) preTest(:,end)];
%myW = prediction(:,[1 7]); % mass volume
%RATIO = mBCC(:,3).*mBCC(:,2).^-1;

%rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23'));
%rmidx = find(strcmp(genoType,'Hp301') | strcmp(genoType,'Il14H') |  strcmp(genoType,'P39') | strcmp(genoType,'bt1/+ in W23') | strcmp(genoType,' bt2/+ in W64A') | strcmp(genoType,' bt1/+ in W64A') | strcmp(genoType,' bt2/+ in W23A'));
%rmidx = prediction(:,1) <100;
%rmidx = RATIO > 700;
%rmidx = RATIO < 50;
%rmidx = RATIO > 1000 | RATIO < 7
%SHOT(rmidx,:) = [];
%subTIP = otipData;
%subTIP(rmidx,:) = [];
%myW(rmidx,:) = [];

[shS shC shU shE shL shERR shLAM] = PCA_FIT_FULL(SHOT,3);
plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.')
hold on
%quiver3(shU(1),shU(2),shU(3),shE(1,1),shE(2,1),shE(3,1),200)
%quiver3(shU(1),shU(2),shU(3),shE(1,1),-shE(2,1),-shE(3,1),200)
%quiver3(shU(1),shU(2),shU(3),shU(1),shU(2),shU(3),2)
quiver3(shU(1),shU(2),shU(3),shU(1),shU(2),shU(3),2)
quiver3(shU(1),shU(2),shU(3),-shU(1),-shU(2),-shU(3),2)
shotVec = shU/norm(shU);
uSHOT = bsxfun(@minus,SHOT,shU);
deltaalong = uSHOT*shotVec';
totalong = uSHOT*shotVec' + norm(shU);
aSHOT = uSHOT - deltaalong*shotVec;
adist = sum(aSHOT.*aSHOT,2).^.5;
figure;
hist(adist,100)
DA = linspace(min(deltaalong),max(deltaalong) - .3*max(deltaalong),100);
da = 20;
h1 = figure;
h2 = figure;
h3 = figure;
plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.');
hold on
vara = [];
ua = [];
Mani = [];
for e =1:numel(DA)    
    fidx = find((deltaalong < DA(e) + da) & (deltaalong > DA(e) - da));
    subA = aSHOT(fidx,:);
    
    subMean = mean(subA,1);
    
    %subA = bsxfun(@minus,subA,subMean);
    
    subA = bsxfun(@plus,subA,shotVec*DA(e)+shU);    
    
    figure(h1);
    plot(subA(:,1),subA(:,2),'.')
    
    
    
    
    [subS subC subU subE subL subERR subLAM] = PCA_FIT_FULL(subA,3);
    theta = linspace(-pi,pi,50);
    dx = 2*(subLAM(1,1).^.5)*cos(theta);
    dy = 2*(subLAM(2,2).^.5)*sin(theta);
    subLAM = diag(diag(subLAM).^-.5);
    subLAM(3,3) = 0;
    
    plot(dx,dy,'r');
    DX = [dx' dy' zeros(size(dx'))];
    Mani(:,:,e) = PCA_BKPROJ(DX,subE,subU);
    
    UMAN(e,:) = subU;
    subE = subE*subLAM;
    EMAN(:,:,e) = subE;
    
    nadist = sum((subC.*subC),2).^.5;
    
    
    
    drawnow
    pause(.05)
    
    vara(e) = std(adist(fidx));
    ua(e) = mean(adist(fidx));
    
    vara(e) = std(nadist);
    ua(e) = mean(nadist);
    
    figure(h2);    
    plot3(SHOT(:,1),SHOT(:,2),SHOT(:,3),'.')
    hold on
    plot3(SHOT(fidx,1),SHOT(fidx,2),SHOT(fidx,3),'r.')
    hold off
    
    figure(h3);    
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
    
    drawnow
end
figure(h3)
Mani = permute(Mani,[3 2 1]);
for e = 1:size(Mani,3)
    plot3(Mani(:,1,e),Mani(:,2,e),Mani(:,3,e),'k')
end
figure;
plot(vara);
figure;
plot(ua);

%% project data to manifold
dU = diff(UMAN,1,1);
dU = cumsum([0;sum(dU.^2,2).^.5]);
for e = 1:size(SHOT,1)
    delta = bsxfun(@minus,UMAN,SHOT(e,:));
    delta = sum(delta.*delta,2);
    [~,midx] = min(delta);
    mySHOT(e,:) = PCA_REPROJ(SHOT(e,:),EMAN(:,:,midx),UMAN(midx,:));
    mySHOT(e,3) = dU(midx);
end
figure;plot3(mySHOT(:,1),mySHOT(:,2),mySHOT(:,3),'.');
%% chekc for  files being there
ext = '_F.tiff';
for e = 1:numel(imageList)
    [I] = readKernelImageFileE(imageList{e},ext);
end
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


for e = 1:10:1000%numel(sMASTERI)
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


for e = 1:10:1000%numel(MASTERI)
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

%% plot lists from NEW-NEW
%mR = mean([RATIO3;RATIO4]);
%mR = nanmean(RATIO3);
close all
[C,ia,ib] = intersect(sideLength(:,2),frontLength(:,2));
pop = pop_front(ib);
%plot(sideLength(ia,1),frontLength(ib,1),'.')
cidx = find(~strcmp(pop,'seedling_phenotyping_widiv'));
widx = find(strcmp(pop,'seedling_phenotyping_widiv'));
subSide = sideLength(ia,:);
subFront = frontLength(ib,:);
%kidx = subFront(find(subSide(:,1) > 700),2);
subFront(cidx,1) = subFront(cidx,1)/nanmean(RATIO);
subSide(cidx,1) = (subSide(cidx,1))/nanmean(RATIO2);
%subFront(widx,1) = subFront(widx,1)/nanmean(mR);
%subSide(widx,1) = subSide(widx,1)/nanmean(mR);
subFront(widx,1) = subFront(widx,1)/nanmean(RATIO4);
subSide(widx,1) = (subSide(widx,1)+21.2)/nanmean(RATIO3);
subFront_Width = frontWidth(ib,:);
subFront_Width(cidx,1) = subFront_Width(cidx,1)/nanmean(RATIO);
%subFront_Width(widx,1) = subFront_Width(widx,1)/nanmean(mR);
subFront_Width(widx,1) = subFront_Width(widx,1)/nanmean(RATIO4);
subSide_Depth = sideDepth(ia,:);
subSide_Depth(cidx,1) = subSide_Depth(cidx,1)/nanmean(RATIO2);
%subSide_Depth(widx,1) = subSide_Depth(widx,1)/mR;
subSide_Depth(widx,1) = subSide_Depth(widx,1)/nanmean(RATIO3);

pops = pop_side(ia);
UQ = unique(pop);
hold on
for u = 1:numel(UQ)
    fidx = strcmp(UQ{u},pop);
    plot3(subSide(fidx,1),subFront(fidx,1),u*ones(1,sum(fidx)),'.')
    hold all
end
plot([0 800],[0 800],'r')
legend(UQ)
oldRes = strcmp(pop,'NAM_parents') | strcmp(pop,'NC-350_RILs');
tofitI = cidx;
p = polyfit(subSide(tofitI,1),subFront(tofitI,1),1);
P = subFront(tofitI,1).*subSide(tofitI,1).^-1;
pv = polyval(p,[0 800]);
plot([0 800],pv,'g')';
%% chekk story re build
cidx = find(~strcmp(pop,'seedling_phenotyping_widiv'));
cidx = 1:numel(pop);
LENGTH = .5*(subFront(cidx,1) + subSide(cidx,1));
WIDTH = subFront_Width(cidx,1);
DEPTH = subSide_Depth(cidx,1);
kid = subFront(cidx,2);


%% get populate for convert to inches for NAM_parents for the front view
oldRes = strcmp(pop_front,'NAM_parents') | strcmp(pop_front,'NC-350_RILs');
oldRes = strcmp(pop_front,'NAM_parents');
kid = frontWidth(oldRes,2);
for e = 1:numel(kid)
    [num2str(frontWidth(e,1)) '-->' num2str(kid(e)) '-->' kernelID.get(num2str(kid(e)))]
end
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
    [num2str(sideDepth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' kernelID.get(num2str(kid(e)))]
end

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
%% get populate for convert to inches for back_room and therefore WIDIV for the side view
oldRes = find(strcmp(pop_side,'seedling_phenotyping_widiv'));
kid = sideDepth(oldRes,2);
for e = 1:numel(kid)
    plt = [num2str(sideDepth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' kernelID.get(num2str(kid(e)))];
    if ~isempty(strfind(plt,'13352'))
        plt
    end
end
%421-->79598-->WISN10_13352*_B_1--seedling_phenotyping_widiv
%554-->79587-->WISN10_13352*_B_4--seedling_phenotyping_widiv
%607-->79589-->WISN10_13352*_B_6--seedling_phenotyping_widiv
%528-->79590-->WISN10_13352*_B_7--seedling_phenotyping_widiv
%% get populate for convert to inches for back_room and therefore WIDIV for the front view
oldRes = find(strcmp(pop_side,'seedling_phenotyping_widiv'));
kid = sideDepth(oldRes,2);
for e = 1:numel(kid)
    plt = [num2str(sideDepth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' kernelID.get(num2str(kid(e)))];
    if ~isempty(strfind(plt,'13352'))
        plt
    end
end
%421-->79598-->WISN10_13352*_B_1--seedling_phenotyping_widiv
%554-->79587-->WISN10_13352*_B_4--seedling_phenotyping_widiv
%607-->79589-->WISN10_13352*_B_6--seedling_phenotyping_widiv
%528-->79590-->WISN10_13352*_B_7--seedling_phenotyping_widiv
%% get populate for convert to inches for back_room and therefore WIDIV for the front view
oldRes = find(strcmp(pop_front,'seedling_phenotyping_widiv'));
kid = frontWidth(oldRes,2);
for e = 1:numel(kid)
    plt = [num2str(frontWidth(oldRes(e),1)) '-->' num2str(kid(e)) '-->' kernelID.get(num2str(kid(e)))];
    if ~isempty(strfind(plt,'13352'))
        plt
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
[c13 c23 V] = impixel(scan2);
WID3 = [c13 c23];
dW3 = diff(WID3,1,1);
dW3 = sum(dW3.*dW3,2).^.5;
dW3 = dW3(1:2:end);
%% make measure by hand and store in hard coded table
TABLE3 = [421;554;607;528];
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
[c14 c24 V] = impixel(scan2);
WID4 = [c14 c24];
dW4 = diff(WID4,1,1);
dW4 = sum(dW4.*dW4,2).^.5;
dW4 = dW4(1:2:end);
%% make measure by hand and store in hard coded table
TABLE4 = [494;538;520;558];
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



%% NEW watch kernel images fly across the sceen from side
ext = '_F.tiff';
PnewF = zeros(numel(f.imageTOP),2);
imageList = f.imageTOP;
clear tmp;
PER = .65;
for e = 1:numel(f.imageTOP)
    try
        try
            tmp = imread([imageList{(e)} ext]);
        catch ME2
            try
                tmp = imread(strrep([imageList{(e)}],'_S2.tiff',ext));
            catch ME2
                try
                    tmp = imread([strrep(imageList{(e)},'gray/07S','gray/07') ext]);
                catch ME2
                    tmp = imread(strrep([strrep(imageList{(e)},'gray/07S','gray/07')],'_S2.tiff',ext));
                end
            end
        end
    catch ME2
    
    end
    
    
    
   tmp(end-100:end,:) = [];
   out = flattenMaskOverlay(tmp,tmp>25);
   imshow(out,[])
   
   drawnow
    e
    numel(imageList)
end
%% load all images



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
[kernelVec genoVec popVec kernelID] = translateWellNames(A);
%save('/mnt/scratch1/phytoM/flashProjects/maize/specKey2.mat','kernelVec','genoVec','popVec','kernelID');
%% try to position to well names
[k g p] = translateWellNames(A);
%% load NAM from database
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...
                'JOIN predictions ' ...
                'on kernels.id = predictions.kernel_id ' ...
                'JOIN nathan_temp ' ...
                'ON kernels.id = nathan_temp.kernel_id ' ...
                'JOIN averageweightspectra_vw ' ...
                'ON kernels.id=averageweightspectra_vw.kernel_id ' ...
                'JOIN kernel_3d ' ...
                'ON kernels.id=kernel_3d.kernel_id ' ...
                'JOIN reporting.root_length_crosstab '...
                'ON kernels.id=reporting.root_length_crosstab.kernel_id '...
                'JOIN reporting.kernel_dims_report_tbl '...
                'on kernels.id=reporting.kernel_dims_report_tbl.kernel_id '...
                'JOIN files ' ...
                'on kernels.img_gray_front_fid = files.id ' ...
                'WHERE population_lines.type = ' '''NAM_parents'' '];
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
%% split data NAM
fidx = [];
fidx(1) = find(strcmp(fieldString,'top_mi_area'));
fidx(2) = find(strcmp(fieldString,'top_major_axis'));
fidx(3) = find(strcmp(fieldString,'top_minor_axis'));

fidx(4) = find(strcmp(fieldString,'front_mi_area'));
fidx(5) = find(strcmp(fieldString,'front_major_axis'));
fidx(6) = find(strcmp(fieldString,'front_minor_axis'));

fidx(7) = find(strcmp(fieldString,'side_mi_area'));
fidx(8) = find(strcmp(fieldString,'side_major_axis'));
fidx(9) = find(strcmp(fieldString,'side_minor_axis'));

OFFSET = 7;
f.kernel_id = cell2mat(results(:,34));
f.tipAngle = cell2mat(results(:,46+OFFSET:46+OFFSET+60));
f.specData = cell2mat(results(:,109+OFFSET:109+OFFSET+781));
f.genoType = results(:,16);
f.predictions = cell2mat(results(:,[35:40 898]));
f.specData = f.specData(:,4:end);
f.shapeData = cell2mat(results(:,fidx));
f.plateName = results(:,6);
f.position = results(:,23);
f.length = cell2mat(results(:,943:943+60));
f.shapeDataL = cell2mat(results(:,1005:1008));
f.xpos = cell2mat(results(:,24));
f.imageTOP = results(:,end);
%% load composition data
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
q = ['SELECT * FROM population_lines JOIN kernel_plates ' ...
    'ON population_lines.id=kernel_plates.population_line_id ' ... 
    'JOIN kernels ' ...
    'ON kernel_plates.id=kernels.plate_id ' ...
    'JOIN predictions ' ...
    'on kernels.id = predictions.kernel_id ' ...
    'JOIN nathan_temp ' ...
    'ON kernels.id = nathan_temp.kernel_id ' ...
    'JOIN averageweightspectra_vw ' ...
    'ON kernels.id=averageweightspectra_vw.kernel_id ' ...
    'JOIN kernel_3d ' ...
    'ON kernels.id=kernel_3d.kernel_id ' ...
    'JOIN reporting.root_length_crosstab '...
    'ON kernels.id=reporting.root_length_crosstab.kernel_id '...
    'JOIN reporting.kernel_dims_report_tbl '...
    'on kernels.id=reporting.kernel_dims_report_tbl.kernel_id '...
    'WHERE population_lines.type = ' '''composition_mutants'' '];
cursor = exec(conn, q);
cursor = fetch(cursor);
fieldString = columnnames(cursor,1);
MUTresults = cursor.Data;
%% split data COMPOSITION MUTS
fidx = [];
fidx(1) = find(strcmp(fieldString,'top_mi_area'));
fidx(2) = find(strcmp(fieldString,'top_major_axis'));
fidx(3) = find(strcmp(fieldString,'top_minor_axis'));

fidx(4) = find(strcmp(fieldString,'front_mi_area'));
fidx(5) = find(strcmp(fieldString,'front_major_axis'));
fidx(6) = find(strcmp(fieldString,'front_minor_axis'));

fidx(7) = find(strcmp(fieldString,'side_mi_area'));
fidx(8) = find(strcmp(fieldString,'side_major_axis'));
fidx(9) = find(strcmp(fieldString,'side_minor_axis'));


OFFSET = 7;

Mf.tipAngle = cell2mat(MUTresults(:,46+OFFSET:46+OFFSET+60));
Mf.specData = cell2mat(MUTresults(:,109+OFFSET:109+OFFSET+781));
Mf.genoType = MUTresults(:,16);
Mf.predictions = cell2mat(MUTresults(:,[35:40 898]));
Mf.specData = Mf.specData(:,4:end);
Mf.shapeData = cell2mat(MUTresults(:,fidx));
Mf.plateName = MUTresults(:,6);
Mf.position = MUTresults(:,23);
Mf.length = cell2mat(MUTresults(:,943:943+60));
Mf.shapeDataL = cell2mat(MUTresults(:,1005:1008));
Mf.xpos = cell2mat(MUTresults(:,24));
Mf.imageTOP = MUTresults(:,end);
%% clean data NAM
GT = 'Hp301';
GT2 = '';
ridx = find(any(abs(diff(f.tipAngle,1,2)) > 10*pi/180,2) | ...
            any(isnan(f.specData),2) | any(isnan(f.tipAngle),2) | ...
            any(abs(f.tipAngle(:,1)) > 45*pi/180,2) | ...
            any(abs(f.tipAngle(:,1) - f.tipAngle(:,end)) < 10*pi/180,2) | ...
            any(f.shapeData(:,4) > 3*10^5,2) | ...
            any(f.shapeData(:,5) > 400,2)  | ...
            any(abs(diff(f.length,1,2)) > 4,2)  | ...
            any(f.shapeDataL(:,end) > .5*10^9) | ...
            f.predictions(:,5) > f.predictions(:,6) | ... 
            isnan(f.xpos) | ...
            any(strcmp(GT2,f.genoType),2) | ...
            any(strcmp(GT,f.genoType),2) ...
            );  
f.tipAngle(ridx,:) = [];
f.kernel_id(ridx,:) = [];

f.specData(ridx,:) = [];
f.genoType(ridx) = [];
f.predictions(ridx,:) = [];
f.shapeData(ridx,:) = [];
f.plateName(ridx,:) = [];
f.position(ridx,:) = [];
f.length(ridx,:) = [];
f.shapeDataL(ridx,:) = [];
f.xpos(ridx,:) = [];
f.imageTOP(ridx) = [];
%% clean data MUTS
ridx = find(any(abs(diff(Mf.tipAngle,1,2)) > 10*pi/180,2) | ...
            any(isnan(Mf.specData),2) | any(isnan(Mf.tipAngle),2) | ...
            any(abs(Mf.tipAngle(:,1)) > 45*pi/180,2) | ...
            any(abs(Mf.tipAngle(:,1) - Mf.tipAngle(:,end)) < 10*pi/180,2) | ...            
            any(Mf.shapeData(:,4) > 3*10^5,2) | ...
            any(Mf.shapeData(:,5) > 400,2)  | ...
            any(abs(diff(Mf.length,1,2)) > 4,2)  | ...
            any(Mf.shapeDataL(:,end) > .5*10^9,2) ...
            ); 
Mf.tipAngle(ridx,:) = [];
Mf.specData(ridx,:) = [];
Mf.genoType(ridx) = [];
Mf.predictions(ridx,:) = [];
Mf.shapeData(ridx,:) = [];
Mf.plateName(ridx,:) = [];
Mf.position(ridx,:) = [];
Mf.length(ridx,:) = [];
Mf.xpos(ridx,:) = [];
Mf.shapeDataL(ridx,:) = [];
Mf.imageTOP(ridx) = [];
%% get labels for both
%LAB.shapeOLD = fieldString(1005:1008);
LAB.shape = {'Height' 'Width' 'Depth' 'Volume'};
LAB.shape = {'Depth' 'Width' 'Height'};
LAB.spec = [fieldString([35:40 898]);'mgOil'];
LAB.spec = {'Percent Protein', 'Percent Oil','Total Density','Material Density','Material Volume','Total Volume','Weight','mg Oil'};
LAB.gravi = {'Swing' 'Final'};
LAB.length = {'Length' 'Final'};
LAB.shapePCA = {'seed set','kernel width','seed size'};
%% get lat para for NAM gravitropism
D = bsxfun(@minus,f.tipAngle,f.tipAngle(:,1));
%D = f.tipAngle;
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL(D,3);
[S1 C1 U1 E1 L1 ERR1 LAM1] = PCA_FIT_FULL(gradient(D),1);
MODEL = C1\wC;
SIM = C1*MODEL;
CS = linspace(min(C1),max(C1),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM,wE,wU);
figure;
subplot(1,2,1);
plot(full_SIM');
DELTA = wC - SIM;
%DELTA = PCA_BKPROJ(DELTA,wE,wU);

[S2 C2 U2 E2 L2 ERR2 LAM2] = PCA_FIT_FULL(DELTA,1);
%check lat para
MODEL2 = C2\wC;
CS = linspace(min(C2),max(C2),10);
SIM2 = CS'*MODEL2;
SIM2 = PCA_BKPROJ(SIM2,wE,wU);
subplot(1,2,2);
plot(SIM2');
f.LAT = [C1 C2];
%% get lat para for MUTS gravitropsim
nD = bsxfun(@minus,Mf.tipAngle,Mf.tipAngle(:,1));
nC1 = PCA_REPROJ(gradient(nD),E1,U1);
wnC = PCA_REPROJ(nD,wE,wU);
SIM = nC1*MODEL;
DELTA2 = wnC - SIM;
nC2 = PCA_REPROJ(DELTA2,E2,U2);
Mf.LAT = [nC1 nC2];
%% get lat para for MUTS length
D = f.length;
[lS lC lU lE lL lERR lLAM] = PCA_FIT_FULL(D,1);
f.GR = lC;
Mf.GR = PCA_REPROJ(Mf.length,lE,lU);
%% add mg oil
[D] = readtext('/mnt/spaldingdata/nate/121213_NIRCalibration_Spectra_NIRTube2_5reps_edited_AVG.csv');
X = cell2mat(D(2:end,26:end));
Y = D(2:end,[15:25]);


toRM = [];
for e = 1:size(Y,1)
    toRM(e) = 0;
    for j = 1:size(Y,2)
        toRM(e) = toRM(e) | isempty(Y{e,j});
    end
end
ridx = find(toRM);
Y(ridx,:) = [];
X(ridx,:) = [];
Y = cell2mat(Y);
dimX = 15;
[caliS caliC caliU caliE caliL caliERR caliLAM] = PCA_FIT_FULL(X,dimX);
dX = mean(X,1);
caliC = bsxfun(@minus,X,mean(X,1));

dY = mean(Y,1);
Y = bsxfun(@minus,Y,dY);
Y = zscore(Y);
perDraw = .01;
clear BETA
for e = 1:size(mean(Y,1),2)
    [Train, Test] = crossvalind('HoldOut', size(caliC,1),perDraw);
    Train = ones(size(Train));
    [XL,YL,XS,YS,BETA(:,e),PCTVAR,MSE,stats] = plsregress(caliC(Train==1,:),Y(Train==1,e),10);
    
    %[A,B] = canoncorr(caliC(Train==1,:),Y(Train==1,e));
    %BETA(2:end,e) = A*inv(B);
    
    yfit = [ones(sum(Train),1) caliC(Train==1,:)]*BETA(:,e);
    ygit = [caliC(Train==1,:)]*BETA(2:end,e);
    
    plot(yfit,Y(Train==1,e),'.');
    
    [R(e),pV] = corr(yfit,Y(Train==1,e));
    title(num2str(R));
    pause(1);
end
%predictVectors = PCA_BKPROJ(BETA(2:end,:)',caliE,caliU);
predictVectors = BETA(1:end,:)';
lessP = predictVectors(:,2:end);

f.predictions = [f.predictions [ones(size(f.specData),1) f.specData]*BETA(:,3)];
%Mf.predictions = [Mf.predictions [ones(size(Mf.specData),1) Mf.specData]*BETA(:,3)];
%% add air space
airSPACE = f.predictions(:,6) - f.predictions(:,5);
f.predictions = [f.predictions airSPACE];
LAB.spec{end+1} = 'airspace';
%% remove air space
f.predictions(:,end) = [];
LAB.spec(end) = [];
%% glue NAM to CCM
M.specData = [f.specData;Mf.specData];
M.predictions = [f.predictions;Mf.predictions];
M.tipAngle = [f.tipAngle;Mf.tipAngle];
M.genoType = [f.genoType;Mf.genoType];
M.shapeData = [f.shapeData;Mf.shapeData];
%M.length = [f.length;Mf.length];
M.LAT = [f.LAT;Mf.LAT];
%M.GR = [f.GR;Mf.GR];
M.plateName = {f.plateName{:} Mf.plateName{:}};
M.position = {f.position{:} Mf.position{:}};
M.shapeDataL = [f.shapeDataL;Mf.shapeDataL];
%% try to position to well names
k = translateWellNames(M);
%% create sense
masselse = f.predictions(:,end-1) - f.predictions(:,1).*f.predictions(:,end-1).^-1 + f.predictions(:,end);
V = f.predictions(:,6).^-1;
M = f.predictions(:,end-1).^-1;
%f.sense = [f.predictions(:,5).*V (f.predictions(:,6) - f.predictions(:,5)).*V f.predictions(:,end).*M (f.predictions(:,1).*f.predictions(:,end-1).^-1).*M masselse.*M];
%f.sense = [f.predictions(:,5).*V (f.predictions(:,6) - f.predictions(:,5)).*V f.predictions(:,2) f.predictions(:,1) masselse.*M];
f.sense = [f.predictions(:,3) f.predictions(:,4) (f.predictions(:,6) - f.predictions(:,5)).*V f.predictions(:,2) f.predictions(:,1) masselse.*M];

%f.sense = [f.predictions(:,end-1).*f.predictions(:,5).^-1 (f.predictions(:,6) - f.predictions(:,5)).*V f.predictions(:,2) f.predictions(:,1) masselse.*M];
%f.sense = [f.predictions(:,4) (f.predictions(:,6) - f.predictions(:,5)).*V f.predictions(:,2) f.predictions(:,1) masselse.*M];

LAB.sense = {'Total Density' 'Material Density' '% Air' '% Oil' '% Protein' '% Else'};
%LAB.sense = {'Material Density' '% Air' '% Oil' '% Protein' '% Else'};
%% ICA
f.ICA = unmix';



%% analysis for FIGURE 3 for paper
T = f;
close all
X = [T.specData];
Y = [T.tipAngle];
%Y = bsxfun(@minus,Y,Y(:,1));
[Sl Xl Ul El Ll ERRl LAMl] = PCA_FIT_FULL(T.length,1);
[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
[S Y Ux Ex L ERR LAM] = PCA_FIT_FULL(Y,3);
Y = [Y Xl];
[Y,mu,sigma] = zscore(Y);
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

for e = 1:size(mU,2)
    figure;
    bar(mB(:,e));
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi));
    title(['Coff Shape:' num2str(e)]);
end

% struct core for Y
[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
    title(['Core Shape:' num2str(e)]);
end

% struct core for X
[Xcore XcoreP] = corr(mU,X);

for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end

for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '-' num2str(mr(e))]);
end
%% curve probe with length
xDIM = 1;
yDIM = 2;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
C = [c1' c2'];
MODEL = C(:,1)\Y;
SIM = C(:,1)*MODEL;
CS = linspace(min(C(:,1)),max(C(:,1)),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM(:,1:3),Ex,Ux);
figure;
subplot(2,2,1);
for e = 1:size(full_SIM,1)
    errorbar(full_SIM(e,:),ER);
    hold all
end




full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM(:,4)*sigma(4),El,Ul);
subplot(2,2,2);
plot(full_SIM');
%%
MODEL = C(:,2)\Y;
SIM = C(:,2)*MODEL;
CS = linspace(min(C(:,2)),max(C(:,2)),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM(:,1:3),Ex,Ux);
subplot(2,2,3);
plot(full_SIM');
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM(:,4),El,Ul);
subplot(2,2,4);
plot(full_SIM');




%% confidance generate ER error
xDIM = 1;
yDIM = 2;
zDIM = 3;
aDIM = 4;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
c3 = [2^-.5 2^-.5]*[mU(:,zDIM) mV(:,zDIM)]';
c4 = [2^-.5 2^-.5]*[mU(:,aDIM) mV(:,aDIM)]';
C = [c1' c2' c3' c4'];
MODEL = C\Y;
SIM = C*MODEL;
MODEL = C\Y;
SIM = C*MODEL;

for tr = 1:size(SIM,1)
    tr_SIM = PCA_BKPROJ(SIM(tr,1:3),Ex,Ux);
    exp_SIM = PCA_BKPROJ(Y(tr,1:3),Ex,Ux);
    differ(tr,:) = abs(tr_SIM - exp_SIM);
    plot(tr_SIM,'r')
    hold on
    plot(exp_SIM,'b')
    hold off
    drawnow
    %pause(.1)
end
ER = mean(differ,1);
%% SUPER BAR good
xDIM = 1;
yDIM = 2;
zDIM = 3;
aDIM = 4;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
c3 = [2^-.5 2^-.5]*[mU(:,zDIM) mV(:,zDIM)]';
c4 = [2^-.5 2^-.5]*[mU(:,aDIM) mV(:,aDIM)]';
C = [c1' c2' c3' c4'];
MODEL = C\Y;
SIM = C*MODEL;


R = corr(C,f.predictions);
for e = 1:2
    figure;
    bar(R(e,:))
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
end

NT = 5;
for i = 1:2
    CS1 = linspace(min(C(:,i)),max(C(:,i)),NT);
    CS2 = linspace(0,0,NT);
    CS3 = linspace(0,0,NT);
    CS4 = linspace(0,0,NT);
    CS = [CS1;CS2;CS3;CS4];
    CS = circshift(CS,[i-1 0]);
    full_SIM = CS'*MODEL;
    full_SIM = bsxfun(@times,full_SIM,sigma);
    full_SIM = PCA_BKPROJ(full_SIM(:,1:3),Ex,Ux);
    
    figure;
    subplot(1,2,1);
    for e = 1:size(full_SIM,1)
        errorbar(180/pi*full_SIM(e,:),180/pi*ER/2);
        hold all
    end
    
    
    csvwrite(['/mnt/spaldingdata/nate/SIM_TA' num2str(i) '.csv'],[full_SIM ;ER]);
   
    %plot(full_SIM');
    full_SIM = CS'*MODEL;
    full_SIM = bsxfun(@times,full_SIM,sigma);
    full_SIM = PCA_BKPROJ(full_SIM(:,4),El,Ul);
    subplot(1,2,2);
    plot(full_SIM');
    %csvwrite('/mnt/spaldingdata/nate/SIM_GR.csv',full_SIM);
end
%% create movies hello fun stuff a must
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
for e = 1:numel(SET)
    [pth nm ext] = fileparts(SET{e}{1});
    fidx = strfind(pth,'/');    
    index{e} = pth(fidx(end)+1:end);
end
%%
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction_13.10.14/loganSPOOL/angle/';
csvFiles = gdig(FilePath,{},{'csv'},1);
DATA = [];
metaData.fileName = {};
metaData.indexName = {};
metaData.index = {};
for e = 1:numel(csvFiles)
    data = csvread(csvFiles{e});
    [pth nm ext] = fileparts(csvFiles{e});
    nm = [nm '_'];
    fidx = strfind(nm,'_');
    fidx(1) = [];
    if size(data,2) == 1
        data = data';
    end
    try
        if size(data,1) == (numel(fidx)-1)
            DATA = [DATA;data];            
            for g = 1:size(data,1)
                metaData.fileName{end+1} = nm;
                metaData.indexName{end+1} = nm(fidx(g)+1:fidx(g+1)-1);
                metaData.index{end+1} = g;
            end
        end        
    end
    e
end
%% search for tip angle movie and seedling
close all
[J qidx] = max(C(:,1));
[J qidx] = max(C(:,2));
%[J qidx] = min(C(:,2));
%[J qidx] = min(C(:,1));

query = f.tipAngle(qidx,:);
an = @(x,y)all(x==y,2);
res = bsxfun(an,DATA,query);
fidx = find(all(res,2));
folderName = metaData.fileName{fidx}(1:end-1);
midx = find(strcmp(folderName,index));
seedlingNumber = metaData.index{midx};
for e = 1:numel(SET{midx})
    filename = SET{midx}{e};
    I = imread(filename);    
    imshow(I)
    title(seedlingNumber)
    drawnow    
end
figure;
plotyy(1:61,query*180/pi,1:61,f.length(qidx,:));
axis([0 61 -30 110]);
%% JUNK SUPER BAR






CS1 = linspace(0,0,10);
CS2 = linspace(min(C(:,2)),max(C(:,2)),10);
CS3 = linspace(0,0,10);
CS4 = linspace(0,0,10);
CS = [CS1;CS2;CS3;CS4];
full_SIM = CS'*MODEL;
full_SIM = bsxfun(@times,full_SIM,sigma);
full_SIM = PCA_BKPROJ(full_SIM(:,1:3),Ex,Ux);
figure;
subplot(1,2,1);
plot(full_SIM');

csvwrite('/mnt/spaldingdata/nate/SIM2.csv',full_SIM);



CS1 = linspace(0,0,10);
CS2 = linspace(0,0,10);
CS3 = linspace(min(C(:,3)),max(C(:,3)),10);
CS = [CS1;CS2;CS3];
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM(:,1:3),Ex,Ux);
figure;
%subplot(2,2,1);
plot(full_SIM');
csvwrite('/mnt/spaldingdata/nate/SIM3.csv',full_SIM);
%% whole structure core
xDIM = 1;
yDIM = 2;
zDIM = 3;
aDIM = 4;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';
c3 = [2^-.5 2^-.5]*[mU(:,zDIM) mV(:,zDIM)]';
c4 = [2^-.5 2^-.5]*[mU(:,aDIM) mV(:,aDIM)]';
C = [c1' c2' c3' c4'];
R = corr(C,f.tipAngle);
Rl = corr(C,f.length);
figure
plot(R(1:2,:)')
figure
plot(Rl(1:2,:)');

%plot(max(gradient(f.tipAngle),[],2),mean(gradient(f.length),2),'.')
%corr(C(:,1),mean(gradient(f.length),2))
%corr(C(:,1),max(gradient(f.tipAngle),[],2))
%% order data by cv1
fidx = find(abs(C(:,2)) < .5);
[JUNK sidx] = sort(C(fidx,1),'descend');
fidx = fidx(sidx);
for e = 1:numel(fidx)
    plot(180/pi*f.tipAngle(fidx(e),:))        
    hold on
    axis([0 61 -20 100])
    pause(.01)
    drawnow
    FA(e) = f.tipAngle(fidx(e),end);
end

%% max variation in CVX


%% create correlation figures for 3 biplot 1-2
xDIM = 1;
yDIM = 2;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';

figure;
hold on
C = [c1' c2'];


LEG = {};
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,1:2);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        %[tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        %plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
        LEG{end+1} = UQ{u};
    end
end


plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');


% try to double check
CHECK = corr(f.predictions,C);
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,1)
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end




% try to double check
GR = mean(gradient(Sl),2);
CHECK = corr([f.LAT GR],C);
L = [LAB.gravi 'GR'];
C_lessP = CHECK;
for e = 1:size(C_lessP,1)    
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
legend(LEG)

axis([-6 6 -6 6])
%% create correlation figures for 3 biplot 2-3

xDIM = 1;
yDIM = 3;
c1 = [2^-.5 2^-.5]*[mU(:,xDIM) mV(:,xDIM)]';
c2 = [2^-.5 2^-.5]*[mU(:,yDIM) mV(:,yDIM)]';


figure;
hold on
C = [c1' c2'];


LEG = {};
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = C(fidx,1:2);
    if size(sub,1) > 10
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        %[tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        %plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
        LEG{end+1} = UQ{u};
    end
end


plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');


% try to double check
CHECK = corr(f.predictions,C);
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,1)
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end


% try to double check
CHECK = corr(f.LAT,C);
L = LAB.gravi;
C_lessP = CHECK;
for e = 1:size(C_lessP,1)    
    P = C_lessP(e,1:2);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
legend(LEG)

axis([-6 6 -6 6])
%% focus simulation on final angle
T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = ([T.tipAngle]);
Y2 = [f.length];
rY = (Y);


[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
%X = [X f.tipAngle(:,1)];
%X = [X X2];
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y,3);
%Y = Y(:,end);

PRE = [];

for u = 1:size(X,1)
    fidx = setdiff(1:size(X,1),u);
    tidx = u;
    
    subX = X(fidx,:);
    subY = Y(fidx,:);
    testX = X(tidx,:);
    testY = Y(tidx,:);
    % network
    net = cascadeforwardnet([10]);
    %net = feedforwardnet([10 3 10]);
    net.trainParam.showWindow = false;
    %net = fitnet([30]);    
    net = train(net,subX',subY');
    netPrediction = net(testX')';
    PRE = [PRE;netPrediction];
    corr(PRE,Y(1:numel(PRE)))
    
end


%% simulation again
T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = ([T.tipAngle]);
%Y = bsxfun(@minus,Y,Y(:,1));
rY = (Y);
%{
% try mutate on the average curve
YU = mean(Y,1);
Z = zeros(size(YU));
for e = 1:numel(YU)
    Z(e,1:e) = YU(1:e);
end
Y = (Z\Y')';
%}

%{
% try spline fit
YC= [ ];
for e = 1:size(Y,1)
    fn{e} = spap2(4,3,1:61,Y(e,:));
    y = fnval(fn{e},1:61);
    %{
    plot(y,'b');
    hold on;
    plot(Y(e,:),'r')
    hold off
    drawnow
    %}
    YC = [YC;fn{e}.coefs];
end
Y = YC;
%}
%{
for e = 1:size(Y,1)
    plot(Z*Y(e,:)','ro');
    hold on
    plot(rY(e,:),'b')
    hold off
    drawnow
    pause(1)
end
%}


Y2 = [f.length];


[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,100);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
%X = [X f.tipAngle(:,end)];
%[X,mu,sigma] = zscore(X);
%X = [X X2];
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
Y = Y(:,29);
Y = mean(gradient(Y),2);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL((Y),1);
COMP_pre = [];
COMP_acc = [];


%{
for u = 1:numel(UQ)
    fidx = strcmp(f.genoType,UQ{u});
    Yo = Y(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Yo,3);
    Y(fidx,:) = bsxfun(@plus,Yo,Uo);
    
    fidx = strcmp(f.genoType,UQ{u});
    Xo = X(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Xo,3);
    X(fidx,:) = bsxfun(@plus,Xo,Uo);
end
%}

for u = 1:numel(UQ)
    fidx = ~strcmp(f.genoType,UQ{u});
    tidx = strcmp(f.genoType,UQ{u});    
    subX = X(fidx,:);
    subY = Y(fidx,:);
    testX = X(tidx,:);
    testY = Y(tidx,:);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(subX,subY);
    % network
    %net = cascadeforwardnet([10]);
    net = feedforwardnet([10]);
    %net = fitnet([30]);    
    net = train(net,subX',subY');
    netPrediction = net(testX')';
    
    % correlation
    testX = bsxfun(@minus,testX,mean(subX,1));
    prediction = (testX*mA)*inv(mB);
    %{
    % correlation + lse
    testX = bsxfun(@minus,testX,mean(subX,1));
    MOD = mV\subY;
    prediction = (testX*mA)*MOD;    
    
    % least squares - from cca
    MODEL = mU\subY;
    netPrediction = (testX*mA)*MODEL;
    % least squares - from raw
    MODEL = subX\subY;
    netPrediction = testX*MODEL;
    %}
    
    %COMP_acc = [COMP_acc;testY];
    %COMP_pre = [COMP_pre;netPrediction];
    %COMP_acc = [COMP_acc;mean(testY,1)];
    %COMP_pre = [COMP_pre;mean(netPrediction,1)];
    
    prediction = PCA_BKPROJ(prediction,Ey,Uy);
    netPrediction = PCA_BKPROJ(netPrediction,Ey,Uy);
    simulation = PCA_BKPROJ(testY,Ey,Uy);
    
    %{
    tmplate = fn{1};
    pre = [];
    npre = [];
    sim = [];
    for e = 1:size(prediction,1)
        tmplate.coefs = prediction(e,:);
        pre(e,:) = fnval(tmplate,1:61);
        tmplate.coefs = netPrediction(e,:);
        npre(e,:) = fnval(tmplate,1:61);        
        tmplate.coefs = simulation(e,:);
        sim(e,:) = fnval(tmplate,1:61);        
    end
    prediction = pre;
    netPrediction = npre;    
    simulation = sim;
    %}
    
    %simulation = simulation*(Z)';
    %netPrediction = netPrediction*(Z)';
    %prediction = prediction*(Z)';
    
    COMP_acc = [COMP_acc;simulation];
    COMP_pre = [COMP_pre;netPrediction];
    
    
    plot(COMP_acc,COMP_pre,'.')
    
    raw = rY(tidx,:);
    %{
    sum((simulation-netPrediction).^2,2).^.5;
    figure;
    for e = 1:size(netPrediction,1)
        plot(simulation(e,:));hold on;plot(netPrediction(e,:),'r');
        hold off
        drawnow
        waitforbuttonpress;
    end
    %}
    
    %{
    errorbar(mean(prediction,1),std(prediction,1,1),'r');
    hold on
    errorbar(mean(netPrediction,1),std(netPrediction,1,1),'b');
    errorbar(mean(simulation,1),std(simulation,1,1),'g');
    errorbar(mean(raw,1),std(raw,1,1),'k');
    drawnow
    hold off
    pause(.5);    
    
    LEG = UQ{u};
    %}
    
end
[COR PVAL] = corr(COMP_acc,COMP_pre);
plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));
plot(COMP_acc,COMP_pre,'.')
%% simulation again
T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = ([T.tipAngle]);

%Y = bsxfun(@minus,Y,Y(:,1));

rY = (Y);
%Y = gradient(Y);
%{
% try mutate on the average curve
YU = mean(Y,1);
Z = zeros(size(YU));
for e = 1:numel(YU)
    Z(e,1:e) = YU(1:e);
end
Y = (Z\Y')';
%}

%{
% try spline fit
YC= [ ];
for e = 1:size(Y,1)
    fn{e} = spap2(4,3,1:61,Y(e,:));
    y = fnval(fn{e},1:61);
    %{
    plot(y,'b');
    hold on;
    plot(Y(e,:),'r')
    hold off
    drawnow
    %}
    YC = [YC;fn{e}.coefs];
end
Y = YC;
%}
%{
for e = 1:size(Y,1)
    plot(Z*Y(e,:)','ro');
    hold on
    plot(rY(e,:),'b')
    hold off
    drawnow
    pause(1)
end
%}


Y2 = [f.length];


[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
%X = [X f.tipAngle(:,1) f.tipAngle(:,29) f.tipAngle(:,end)];
%X = [f.tipAngle(:,1:2)  mean(gradient(f.tipAngle),2)];
%X = [f.tipAngle(:,1) f.tipAngle(:,29) f.tipAngle(:,end)];
%X = [X2];
[X,mu,sigma] = zscore(X);
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL((Y),3);
COMP_pre = [];
COMP_acc = [];


%{
for u = 1:numel(UQ)
    fidx = strcmp(f.genoType,UQ{u});
    Yo = Y(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Yo,3);
    Y(fidx,:) = bsxfun(@plus,Yo,Uo);
    
    fidx = strcmp(f.genoType,UQ{u});
    Xo = X(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Xo,3);
    X(fidx,:) = bsxfun(@plus,Xo,Uo);
end
%}

for u = 1:numel(UQ)
    fidx = ~strcmp(f.genoType,UQ{u});
    tidx = strcmp(f.genoType,UQ{u});    
    subX = X(fidx,:);
    subY = Y(fidx,:);
    testX = X(tidx,:);
    testY = Y(tidx,:);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(subX,subY);
    % network
    %net = cascadeforwardnet([10]);
    net = feedforwardnet([3]);
    %net = fitnet([30]);    
    %net = train(net,subX',subY');
    %netPrediction = net(testX')';
    
    % correlation
    testX = bsxfun(@minus,testX,mean(subX,1));
    prediction = (testX*mA)*inv(mB);
    %{
    % correlation + lse
    testX = bsxfun(@minus,testX,mean(subX,1));
    MOD = mV\subY;
    prediction = (testX*mA)*MOD;    
    
    % least squares - from cca
    MODEL = mU\subY;
    netPrediction = (testX*mA)*MODEL;
    % least squares - from raw
    MODEL = subX\subY;
    netPrediction = testX*MODEL;
    %}
    
    %COMP_acc = [COMP_acc;testY];
    %COMP_pre = [COMP_pre;netPrediction];
    %COMP_acc = [COMP_acc;mean(testY,1)];
    %COMP_pre = [COMP_pre;mean(netPrediction,1)];
    
    prediction = PCA_BKPROJ(prediction,Ey,Uy);
    %netPrediction = PCA_BKPROJ(netPrediction,Ey,Uy);
    simulation = PCA_BKPROJ(testY,Ey,Uy);
    
    %simulation = cumsum(simulation,2);
    %netPrediction = cumsum(netPrediction,2);
    %prediction = cumsum(prediction,2);
    
    %{
    tmplate = fn{1};
    pre = [];
    npre = [];
    sim = [];
    for e = 1:size(prediction,1)
        tmplate.coefs = prediction(e,:);
        pre(e,:) = fnval(tmplate,1:61);
        tmplate.coefs = netPrediction(e,:);
        npre(e,:) = fnval(tmplate,1:61);        
        tmplate.coefs = simulation(e,:);
        sim(e,:) = fnval(tmplate,1:61);        
    end
    prediction = pre;
    netPrediction = npre;    
    simulation = sim;
    %}
    
    %simulation = simulation*(Z)';
    %netPrediction = netPrediction*(Z)';
    %prediction = prediction*(Z)';
    
    COMP_acc = [COMP_acc;mean(simulation,1)];
    COMP_pre = [COMP_pre;mean(prediction,1)];
    
    raw = rY(tidx,:);
    %{
    sum((simulation-netPrediction).^2,2).^.5;
    figure;
    for e = 1:size(netPrediction,1)
        plot(simulation(e,:));hold on;plot(netPrediction(e,:),'r');
        hold off
        drawnow
        waitforbuttonpress;
    end
    %}
    
    
    errorbar(mean(prediction,1),std(prediction,1,1),'r');
    hold on
    errorbar(mean(netPrediction,1),std(netPrediction,1,1),'b');
    errorbar(mean(simulation,1),std(simulation,1,1),'g');
    errorbar(mean(raw,1),std(raw,1,1),'k');
    drawnow
    hold off
    pause(.5);    
    
    LEG = UQ{u};
    
    
end
[COR PVAL] = corr(COMP_acc,COMP_pre);
%sCOR = diag(COR);
plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));
%plot(COMP_acc(:,29),COMP_pre(:,29),'.')
%% simulation again2

T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = [T.tipAngle];
Y2 = [f.length];
%Y = bsxfun(@minus,Y,Y(:,1));
rY = (Y);
[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,25);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
X = [X X2];
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y,3);
COMP_pre = [];
COMP_acc = [];

%{
for u = 1:numel(UQ)
    fidx = strcmp(f.genoType,UQ{u});
    Yo = Y(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Yo,3);
    Y(fidx,:) = bsxfun(@plus,Yo,Uo);
    
    fidx = strcmp(f.genoType,UQ{u});
    Xo = X(fidx,:);
    [So Yo Uo Eo Lo ERRo LAMo] = PCA_FIT_FULL(Xo,3);
    X(fidx,:) = bsxfun(@plus,Xo,Uo);
end
%}

[EXP PVAL] = probeTopology(@(X)cascadeforwardnet(X),10,3,3,X,Y,T.genoType,Ey,Uy);
[COR PVAL] = corr(COMP_acc,COMP_pre);
plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL))
%% simulation again hold single out

T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = gradient([T.tipAngle]);
Y2 = [f.length];
%Y = bsxfun(@minus,Y,Y(:,1));





%{
UQ = unique(f.genoType)
for u = 1:numel(UQ)
    fidx = strcmp(UQ{u},f.genoType);
    Y(fidx,:) = bsxfun(@minus,Y(fidx,:),mean(Y(fidx,:)));
    X(fidx,:) = bsxfun(@minus,X(fidx,:),mean(X(fidx,:)));    
end
rY = (Y);
%}


[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
%X = [X f.tipAngle(:,1)];
%X = [X X2];
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y,3);





COMP_pre = [];
COMP_pre_net = [];
COMP_acc = [];

for u = 1:size(T.tipAngle,1)
    fidx = setdiff(1:size(T.tipAngle,1),u);
    tidx = u;
    
    subX = X(fidx,:);
    subY = Y(fidx,:);
    testX = X(tidx,:);
    testY = Y(tidx,:);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(subX,subY);
    % network
    %net = cascadeforwardnet([10]);
    %net = feedforwardnet([10 3 10]);
    %net.trainParam.showWindow = false;
    %net = fitnet([30]);    
    %net = train(net,subX',subY');
    %netPrediction = net(testX')';
    
    % correlation
    testX = bsxfun(@minus,testX,mean(subX,1));
    prediction = (testX*mA)*inv(mB);
    %prediction = (testX*mA);
    %{
    % correlation + lse
    testX = bsxfun(@minus,testX,mean(subX,1));
    MOD = mV\subY;
    prediction = (testX*mA)*MOD;    
    
    % least squares - from cca
    MODEL = mU\subY;
    netPrediction = (testX*mA)*MODEL;
    % least squares - from raw
    MODEL = subX\subY;
    netPrediction = testX*MODEL;
    %}
    
    %COMP_acc = [COMP_acc;testY];
    %COMP_pre = [COMP_pre;netPrediction];
    %COMP_acc = [COMP_acc;mean(testY,1)];
    %COMP_pre = [COMP_pre;mean(netPrediction,1)];
    
    prediction = PCA_BKPROJ(prediction,Ey,Uy);
    %netPrediction = PCA_BKPROJ(netPrediction,Ey,Uy);
    simulation = PCA_BKPROJ(testY,Ey,Uy);
    
    %{
    tmplate = fn{1};
    pre = [];
    npre = [];
    sim = [];
    for e = 1:size(prediction,1)
        tmplate.coefs = prediction(e,:);
        pre(e,:) = fnval(tmplate,1:61);
        tmplate.coefs = netPrediction(e,:);
        npre(e,:) = fnval(tmplate,1:61);        
        tmplate.coefs = simulation(e,:);
        sim(e,:) = fnval(tmplate,1:61);        
    end
    prediction = pre;
    netPrediction = npre;    
    simulation = sim;
    %}
    
    %simulation = simulation*(Z)';
    %netPrediction = netPrediction*(Z)';
    %prediction = prediction*(Z)';
    
    COMP_acc = [COMP_acc;simulation];
    COMP_pre = [COMP_pre;prediction];
    %COMP_pre_net = [COMP_pre;netPrediction];
    
    raw = rY(tidx,:);
    %{
    sum((simulation-netPrediction).^2,2).^.5;
    figure;
    for e = 1:size(netPrediction,1)
        plot(simulation(e,:));hold on;plot(netPrediction(e,:),'r');
        hold off
        drawnow
        waitforbuttonpress;
    end
    %}
    
    %{
    errorbar(180/pi*mean(prediction,1),std(prediction,1,1),'r');
    hold on
    %errorbar(180/pi*mean(netPrediction,1),std(netPrediction,1,1),'b');
    errorbar(180/pi*mean(simulation,1),std(simulation,1,1),'g');
    errorbar(180/pi*mean(raw,1),std(raw,1,1),'k');
    axis([0 61 -20 120]);
    %axis([0 61 0 800]);
    drawnow
    hold off
    pause(.1);    
    %}
    u
    %LEG = UQ{u};
    
    
end

[COR PVAL] = corr(COMP_acc,COMP_pre);
%sCOR = diag(COR);
figure;plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));

[COR PVAL] = corr(COMP_acc,repmat(mean(rY,1),[size(rY,1) 1]));
%sCOR = diag(COR);
figure;plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));
TIPA = cumsum(COMP_pre,2);
L = corr(TIPA,T.tipAngle);
figure;plot(diag(L))
%% simulation again hold single out and clamp at genotypes means

T = f;
UQ = unique(T.genoType);
X = [T.specData];
X2 = [T.shapeDataL];
Y = ([T.tipAngle]);
Y2 = [f.length];





rY = (Y);


[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
[S X2 Ux2 Ex2 L ERR LAM] = PCA_FIT_FULL(X2,3);
%X = [X f.tipAngle(:,1)];
%X = [X X2];
%[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y2,1);
[S Y Uy Ey L ERR LAM] = PCA_FIT_FULL(Y,3);




UQ = unique(f.genoType);
clear uX uY;
for u = 1:numel(UQ)
    fidx = strcmp(UQ{u},f.genoType);
    uX(u,:) = mean(X(fidx,:));
    uY(u,:) = mean(Y(fidx,:));    
end
X = uX;
Y = uY;




COMP_pre = [];
COMP_pre_net = [];
COMP_acc = [];

for u = 1:size(X,1)
    fidx = setdiff(1:size(X,1),u);
    tidx = u;
    
    subX = X(fidx,:);
    subY = Y(fidx,:);
    testX = X(tidx,:);
    testY = Y(tidx,:);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(subX,subY);
    % network
    %net = cascadeforwardnet([10]);
    net = feedforwardnet([10 3 10]);
    net.trainParam.showWindow = false;
    %net = fitnet([30]);    
    net = train(net,subX',subY');
    netPrediction = net(testX')';
    
    % correlation
    testX = bsxfun(@minus,testX,mean(subX,1));
    prediction = (testX*mA)*inv(mB);
    %prediction = (testX*mA);
    %{
    % correlation + lse
    testX = bsxfun(@minus,testX,mean(subX,1));
    MOD = mV\subY;
    prediction = (testX*mA)*MOD;    
    
    % least squares - from cca
    MODEL = mU\subY;
    netPrediction = (testX*mA)*MODEL;
    % least squares - from raw
    MODEL = subX\subY;
    netPrediction = testX*MODEL;
    %}
    
    %COMP_acc = [COMP_acc;testY];
    %COMP_pre = [COMP_pre;netPrediction];
    %COMP_acc = [COMP_acc;mean(testY,1)];
    %COMP_pre = [COMP_pre;mean(netPrediction,1)];
    
    prediction = PCA_BKPROJ(prediction,Ey,Uy);
    netPrediction = PCA_BKPROJ(netPrediction,Ey,Uy);
    simulation = PCA_BKPROJ(testY,Ey,Uy);
    
    %{
    tmplate = fn{1};
    pre = [];
    npre = [];
    sim = [];
    for e = 1:size(prediction,1)
        tmplate.coefs = prediction(e,:);
        pre(e,:) = fnval(tmplate,1:61);
        tmplate.coefs = netPrediction(e,:);
        npre(e,:) = fnval(tmplate,1:61);        
        tmplate.coefs = simulation(e,:);
        sim(e,:) = fnval(tmplate,1:61);        
    end
    prediction = pre;
    netPrediction = npre;    
    simulation = sim;
    %}
    
    %simulation = simulation*(Z)';
    %netPrediction = netPrediction*(Z)';
    %prediction = prediction*(Z)';
    
    COMP_acc = [COMP_acc;simulation];
    COMP_pre = [COMP_pre;prediction];
    %COMP_pre_net = [COMP_pre;netPrediction];
    
    raw = rY(tidx,:);
    %{
    sum((simulation-netPrediction).^2,2).^.5;
    figure;
    for e = 1:size(netPrediction,1)
        plot(simulation(e,:));hold on;plot(netPrediction(e,:),'r');
        hold off
        drawnow
        waitforbuttonpress;
    end
    %}
    
    
    errorbar(180/pi*mean(prediction,1),std(prediction,1,1),'r');
    hold on
    errorbar(180/pi*mean(netPrediction,1),std(netPrediction,1,1),'b');
    errorbar(180/pi*mean(simulation,1),std(simulation,1,1),'g');
    errorbar(180/pi*mean(raw,1),std(raw,1,1),'k');
    axis([0 61 -20 120]);
    %axis([0 61 0 800]);
    drawnow
    hold off
    pause(.1);    
    
    u
    %LEG = UQ{u};
    
    
end

[COR PVAL] = corr(COMP_acc,COMP_pre);
%sCOR = diag(COR);
figure;plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));

[COR PVAL] = corr(COMP_acc,repmat(mean(rY,1),[size(rY,1) 1]));
%sCOR = diag(COR);
figure;plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));
%%
UQ = unique(f.genoType)
figure;
clear LEG
for u = 1:numel(UQ)
    fidx = strcmp(UQ{u},f.genoType);
    [COR PVAL] = corr(COMP_acc(fidx,:),COMP_pre(fidx,:));
    %plotyy(1:size(COR,1),diag(COR),1:size(COR,1),diag(PVAL));
    plot(1:size(COR,1),diag(COR));
    hold all
    LEG{u} = UQ{u};
end
legend(LEG)

%% rank cv1 and 2
T = f;
fidx = find(abs(C(:,2)) < .5);
[JUNK sidx] = sort(C(fidx,1),'descend');
fidx = fidx(sidx);
fidx = fidx(1:6:end);
figure;
[SIM Cx Ux Ex L ERR LAM] = PCA_FIT_FULL(gradient(T.tipAngle),3);
for i = 1:numel(fidx)
    %ta(i,:) = T.LAT(fidx(i),:);
    
end
plot(ta);
legend(LAB.gravi)
%% match A with new shape from above ZIP WHOLE
nA = A;
nA.SHAPE = [];
toR = [];
for e = 1:numel(nA.kernel_id)
    k = nA.kernel_id{e};
    fidx = find(kid == k);
    if isempty(fidx)
        toR = [toR e];
    else
        nA.SHAPE = [nA.SHAPE;[DEPTH(fidx) WIDTH(fidx) LENGTH(fidx)]];
    end
end
nA.specData(toR,:) = [];
preTest2(toR,:) = []; % obtained from above
nA.plateName(toR) = [];
nA.kernel_id(toR) = [];
%% A) new spec to shape cca
T = nA;%f;
%pidx = find(~strcmp(T.genoType,'W22^ACR') & ~strcmp(T.genoType,'P39') & ~strcmp(T.genoType,'Hp301'));
%pidx = find(~strcmp(T.genoType,'Hp301'));
%pidx = readtext('/home/nate/Downloads/out.csv');
%pidx = cell2mat(pidx(2:end,2));
%[~,pidx] = setdiff(T.kernel_id,pidx);
%pidx = 1:size(preTest,1);
close all
%X = [T.specData];
%X = preTest;
%X = T.specData;
%subIDX = [1 4 5 7 9];% 12];
%X = X(:,subIDX);
X = T.specData(:,4:end);
%X = preTest2(:,subIDX);
%X = preTest;
%Y = [T.shapeDataL(:,1:3)];
Y = [T.SHAPE];
%Y = unmix';
%X = X(pidx,:);
%Y = Y(pidx,:);
nY = Y;
kpidx = (Y(:,2) < 900) & (Y(:,3) < 800);
X = X(kpidx,:);
Y = Y(kpidx,:);



%X = readtext('/home/nate/Downloads/out.NIR.csv');
%Y = readtext('/home/nate/Downloads/out.shape.csv');
%Y = cell2mat(Y(2:end,2:end));


%X = exp(X).^-1;
%X = zscore(X);
%Y = zscore(Y);

[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,15);

%X = zscore(X);
%Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    %set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

for e = 1:size(mV,2)
    figure;
    bar(mB(:,e));
    %set(gca,'XTickLabel',LAB.ICA,'XTick',1:numel(LAB.ICA));
    title(['Coff Shape:' num2str(e)]);
end


[Ycore YcoreP] = corr(mV,Y);
%[Ycore YcoreP] = corr(mU,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.shape(1:3),'XTick',1:numel(LAB.shape(1:3)))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
[Xcore XcoreP] = corr(mU,preTest2(kpidx,subIDX));
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    %set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    set(gca,'XTickLabel',NMS(subIDX),'XTick',1:numel(subIDX))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(corr(mU(:,e),mV(:,e)))]);
end
%% genotype hold out X number of groups YES THIS ONE
close all
T = nA;%f;

%X = preTest2(:,subIDX);
%X = preTest;
%Y = [T.shapeDataL(:,1:3)];

%Y = grData(:,end);


%[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,3);


%Y = bsxfun(@minus,Y,mean(Y,1));
%X = bsxfun(@minus,X,mean(X,1));
%{
for e = 1:size(Y,1)
    Y(e,:) = Y(e,:) * norm(Y(e,:))^-1;
end
for e = 1:size(X,1)
    X(e,:) = X(e,:) * norm(X(e,:))^-1;
end
%}
UQ = unique(T.plateName);

fG = [];
fG2 = [];
fG3 = [];

MASFG = [];
MASFG2 = [];
MAD = [];
UMT = [];
UMB = [];
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
MC = [];
for factors = 3:600
    Y = [T.SHAPE];
    X = T.specData(:,4:end);
    [S X U Ex L ERR LAM] = PCA_FIT_FULL(X,factors);
    for pull = 1:50
        fG = [];
        fG2 = [];
        fG3 = [];
        R = randperm(numel(UQ));
        R = R(1:50);
        sidx0 = [];
        sidx1 = [];
        for e = 1:numel(R)
            %sidx0 = [sidx0 find(~strcmp(pl,UQ{R(e)}))];
            sidx1 = [sidx1;find(strcmp(T.plateName,UQ{R(e)}))];
        end
        sidx0 = setdiff(1:numel(T.plateName),sidx1);
        trainX = X(sidx0,:);
        trainY = Y(sidx0,:);
        testX = X(sidx1,:);        
        %testX = mean(testX,1);
        testY = Y(sidx1,:);
        %testY = mean(testY,1);

        uSx = mean(trainX);
        uSy = mean(trainY);
    
        [A,B,r,trainX,trainY,stats] = canoncorr(trainX,trainY);
        testX = bsxfun(@minus,testX,uSx);
        testY = bsxfun(@minus,testY,uSy);
        testX = testX*A;
        testY = testY*B;
        [tC] = corr(testX,testY);
        if any(isnan(tC(:)))
            report = 1;
            holVx = testX;
            holVy = testY;
            
            break;
        end
        MC(pull,:,factors-2) = diag(tC);
        MCMOD(pull,:,factors-2) = r';
    end
    factors;
    uMC = mean(MC,1);
    uMCMOD = mean(MCMOD,1);
    plot(squeeze(uMC)')
    hold on
    plot(squeeze(uMCMOD)','--')
    hold off
    drawnow
end
%%
WV = csvread('~/T.csv');
newMM = bsxfun(@minus,f.specData(:,4:end),WV(:,1)');
newMM = newMM*WV(:,2:end);
%% Z) added new for check with new data
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle,5);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData,10);
close all
X = sC;
Y = tC;
%X = zscore(X);
%Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    %set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Shape:' num2str(e)]);
end


[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mV,2)
    figure;
    bar(Ycore(e,:));
    %set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi));
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mU,2)
    figure;
    bar(Xcore(e,:));
    %set(gca,'XTickLabel',[LAB.spec],'XTick',1:numel([LAB.spec]))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mU,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% A) no hold out for [raw_spec] to [ICA_seedshape]
[unmix, A, W] = fastica(T.shapeDataL');
%%
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];


X = zscore(X);
Y = zscore(Y);

[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,8);
%[S Y U Ey L ERR LAM] = PCA_FIT_FULL(Y,3);

%{
for e = 1:size(Ex,2)
    figure;
    bar(Ex(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

for e = 1:size(Ey,2)
    figure;
    bar(Ey(:,e));
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Spec:' num2str(e)]);
end
%}

X = zscore(X);
Y = zscore(Y);


[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);


for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

for e = 1:size(mV,2)
    figure;
    bar(mB(:,e));
    set(gca,'XTickLabel',LAB.ICA,'XTick',1:numel(LAB.ICA));
    title(['Coff Shape:' num2str(e)]);
end


[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.shapePCA,'XTick',1:numel(LAB.shapePCA))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end




for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end



figure;
[JUNK meanidx] = min(abs(mU(:,2) - mean(mU(:,2))));
meanImage = imread([T.imageTOP{meanidx} '_T.tiff']);
[JUNK minidx] = min(mU(:,2));
minImage = imread(strrep([T.imageTOP{minidx}],'_F','_T'));
[JUNK maxidx] = max(mU(:,2));
maxImage = imread(strrep([T.imageTOP{maxidx}],'_F','_T'));
imshow([minImage meanImage maxImage],[])

[JUNK meanidx] = min(abs(mU(:,2) - mean(mU(:,2))));
meanImage = imread([T.imageTOP{meanidx} '_F.tiff']);
[JUNK minidx] = min(mU(:,2));
minImage = imread(strrep([T.imageTOP{minidx}],'_F','_F'));
[JUNK maxidx] = max(mU(:,2));
maxImage = imread(strrep([T.imageTOP{maxidx}],'_F','_F'));
imshow([minImage meanImage maxImage],[])
%% A) no hold out for [raw_spec] to [gravi]
T = f;
close all
X = [T.specData];
Y = [T.tipAngle];

[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,30);
[S Y U Ex L ERR LAM] = PCA_FIT_FULL(Y,5);

X = zscore(X);
Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

UQ = unique(T.genoType);
%% A) no hold out for [composition_sense] to [shape]
T = f;
close all
X = [T.sense];
Y = [T.shapeDataL(:,1:3)];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

[wx,wy,pX,pY,D] = myCCA(X,Y,3);
percentX = std(mU,1,1)/sum(std(X,1,1));
percentY = std(mV,1,1)/sum(std(Y,1,1));
XtoY = mA*(mV\Y);
preY = X*XtoY;
perY = var(preY,1,1)/sum(var(Y,1,1));
[S C U E L ERR LAM] = PCA_FIT_FULL(X,size(X,2));

[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.shape(1:3),'XTick',1:numel(LAB.shape(1:3)))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',LAB.sense,'XTick',1:numel(LAB.sense))
    title(['Core Spec:' num2str(e)]);
end

for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% try gravi regression along packing
D = T.length;
e=2;
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL(D,3);
MODEL = mU(:,e)\wC;
CS = linspace(min(mU(:,e)),max(mU(:,e)),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM,wE,wU);
figure;
plot(full_SIM');
(full_SIM(1,end) - full_SIM(end,end))/full_SIM(end,end)

MODEL = mV(:,e)\wC;
CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
full_SIM = CS'*MODEL;
full_SIM = PCA_BKPROJ(full_SIM,wE,wU);
figure;
plot(full_SIM');

(full_SIM(1,end) - full_SIM(end,end))/full_SIM(end,end)
%% B) no hold out for [composition;shape] -> gravi
%B =f
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:4)];
X = [X Y];
Z = [T.LAT];
Y = Z;
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',[LAB.spec;LAB.shape(1:4)'],'XTick',1:numel([LAB.spec;LAB.shape(1:4)']))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end

UQ = unique(T.genoType);
nV = [];
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    nV(u) = mean(mU(fidx,2));
end
[jV sidx] = sort(nV);
sUQ = UQ(sidx);
figure;
bar(jV)
set(gca,'XTickLabel',[sUQ],'XTick',1:numel(sUQ))



figure;
[JUNK meanidx] = min(abs(mU(:,2) - mean(mU(:,2))));
meanImage = imread([T.imageTOP{meanidx} '_T.tiff']);
[JUNK minidx] = min(mU(:,2));
minImage = imread(strrep([T.imageTOP{minidx}],'_F','_T'));
[JUNK maxidx] = max(mU(:,2));
maxImage = imread(strrep([T.imageTOP{maxidx}],'_F','_T'));
imshow([minImage meanImage maxImage],[])
figure;
[JUNK meanidx] = min(abs(mU(:,2) - mean(mU(:,2))));
meanImage = imread([T.imageTOP{meanidx} '_F.tiff']);
[JUNK minidx] = min(mU(:,2));
minImage = imread(strrep([T.imageTOP{minidx}],'_F','_F'));
[JUNK maxidx] = max(mU(:,2));
maxImage = imread(strrep([T.imageTOP{maxidx}],'_F','_F'));
imshow([minImage meanImage maxImage],[])
%% C) no hold out for [shape] -> [gravi]
T = f;
close all
X = [T.shapeDataL(:,1:4)];
Y = [T.LAT];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',[LAB.shape(1:4)],'XTick',1:numel([LAB.shape(1:4)]))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% C) no hold out for [composition] -> gravi
T = f;
close all
X = [T.predictions];
Y = [T.LAT];
%Y = unmix';
%[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,size(X,2));
%[S Y U Ey L ERR LAM] = PCA_FIT_FULL(Y,size(Y,2));
X = zscore(X);
Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
%nX = bsxfun(@minus,X,mean(X));
%nU = nX*mA(:,1);
for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Shape:' num2str(e)]);
end


[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mV,2)
    figure;
    bar(Ycore(e,:));
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi));
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mU,2)
    figure;
    bar(Xcore(e,:));
    set(gca,'XTickLabel',[LAB.spec],'XTick',1:numel([LAB.spec]))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mU,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% C) no hold out for [composition_sense] to [gravi]
T = f;
close all
X = [T.sense];
Y = [T.LAT];
X = zscore(X);
Y = zscore(Y);

[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

for e = 1:size(mV,2)
    figure;
    bar(mA(:,e));
    set(gca,'XTickLabel',LAB.sense,'XTick',1:numel(LAB.sense));
    title(['Coff Shape:' num2str(e)]);
end

[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi));
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',LAB.sense,'XTick',1:numel(LAB.sense))
    title(['Core Spec:' num2str(e)]);
end

for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% C) no hold out for [composition] -> [length gravi]
T = f;
close all
X = [T.predictions];
Y = [T.LAT T.GR];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

%{
[XL,YL,mU,mV,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,2);
mA = X\XS;
mB = Y\YS;
%}


[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mV,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',[LAB.gravi LAB.length],'XTick',1:numel([LAB.gravi LAB.length]))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mU,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',[LAB.spec],'XTick',1:numel([LAB.spec]))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mU,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% C) no hold out for [composition_sense] -> [length gravi]
T = f;
close all
X = [T.sense];
Y = [T.LAT T.GR];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);

%{
[XL,YL,mU,mV,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,2);
mA = X\XS;
mB = Y\YS;
%}


[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mV,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',[LAB.gravi LAB.length],'XTick',1:numel([LAB.gravi LAB.length]))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mU,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',[LAB.sense],'XTick',1:numel([LAB.sense]))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mU,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end
%% B) no hold out for [composition;shape] -> [gravi length]
%B =f
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];
X = [X Y];
Z = [T.LAT T.GR];
Y = Z;
X = zscore(X);
Y = zscore(Y);
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
mA(:,1) = -mA(:,1);
mU = bsxfun(@minus,X,mean(X,1))*mA;
[Ycore YcoreP] = corr(mV,Y);
for e = 1:size(mB,2)
    figure;
    bar(Ycore(e,:))
    set(gca,'XTickLabel',[LAB.gravi LAB.length],'XTick',1:numel([LAB.gravi LAB.length]))
    title(['Core Shape:' num2str(e)]);
end

[Xcore XcoreP] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(Xcore(e,:))
    set(gca,'XTickLabel',[LAB.spec';LAB.shape(1:3)'],'XTick',1:numel([LAB.spec';LAB.shape(1:3)']))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e) '--' num2str(mr(e))]);
end

UQ = unique(T.genoType);
nV = [];
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    nV(u) = mean(mU(fidx,2));
end
[jV sidx] = sort(nV);
sUQ = UQ(sidx);
figure;
bar(jV)
set(gca,'XTickLabel',[sUQ],'XTick',1:numel(sUQ))
%% D) HOLD OUT for [composition] to [shape]
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];
for e = 1:min(size(X,2),size(Y,2))
    fig_ycore(e) = figure;
    fig_xcore(e) = figure;
    over_all(e) = figure;
end

%X = zscore(X);
%Y = zscore(Y);
perDraw = .5;
Ycore = [];
Xcore = [];
mr = [];
for tr = 1:100
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(X,1),perDraw);
    [mA,mB,mr(tr,:),mU,mV,mstats] = canoncorr(X(TrainMaster,:),Y(TrainMaster,:));
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(TrainMaster,:),Y(TrainMaster,:));
    %rho = corr(XS,X(TrainMaster,:));
    
    [Ycore(:,:,tr) YcoreP] = corr(mV,Y(TrainMaster,:));
    if tr > 1
        uYcore = mean(Ycore(:,:,1:end-1),3);
        for can_var = 1:size(uYcore,1)
            if sign(uYcore(can_var,:)*Ycore(can_var,:,tr)') == -1
                Ycore(can_var,:,tr) = -Ycore(can_var,:,tr);
            end
        end
    end
    uYcore = mean(Ycore,3);
    sdYcore = std(Ycore,1,3);
    for e = 1:size(mB,2)
        figure(fig_ycore(e))        
        bar(uYcore(e,:));
        hold on
        h = errorbar(uYcore(e,:),sdYcore(e,:));        
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',LAB.shape(1:3),'XTick',1:numel(LAB.shape(1:3)))
        title(['Canonical Loadings:' num2str(e)]);
        ylabel('Correlation');
        hold off
        drawnow
    end

    [Xcore(:,:,tr) XcoreP] = corr(mU,X(TrainMaster,:));
    if tr > 1
        uXcore = mean(Xcore(:,:,1:end-1),3);
        for can_var = 1:size(uXcore,1)
            if sign(uXcore(can_var,:)*Xcore(can_var,:,tr)') == -1
                Xcore(can_var,:,tr) = -Xcore(can_var,:,tr);
            end
        end
    end
    uXcore = mean(Xcore,3);
    sdXcore = std(Xcore,1,3);
    for e = 1:size(mA,2)
        figure(fig_xcore(e));
        bar(uXcore(e,:));
        hold on
        h = errorbar(uXcore(e,:),sdXcore(e,:));
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
        ylabel('Correlation');
        title(['Canonical Loadings:' num2str(e)]);
        hold off
    end


    for e = 1:size(mA,2)
        figure(over_all(e))
        plot(mU(:,e),mV(:,e),'.')
        hold on
        title(['Overall:' num2str(e)]);
    end

end
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_shape_average_coreX.csv',mean(Xcore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_shape_std_coreX.csv',std(Xcore,0,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_shape_average_coreY.csv',mean(Ycore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_shape_std_coreY.csv',std(Ycore,0,3));
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/lat_core_struct_predictions.csv',s_struct_C);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var1.csv',[mU(:,1),mV(:,1)]);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var2.csv',[mU(:,2),mV(:,2)]);
%% E) HOLD OUT for [composition;shape] -> gravi
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];
X = [X Y];
Z = [T.LAT];
Y = Z;
for e = 1:min(size(X,2),size(Y,2))
    fig_ycore(e) = figure;
    fig_xcore(e) = figure;
    over_all(e) = figure;
end


perDraw = .5;
Ycore = [];
Xcore = [];
for tr = 1:100
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(X,1),perDraw);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(X(TrainMaster,:),Y(TrainMaster,:));
    [Ycore(:,:,tr) YcoreP] = corr(mV,Y(TrainMaster,:));
    if tr > 1
        uYcore = mean(Ycore(:,:,1:end-1),3);
        for can_var = 1:size(uYcore,1)
            if sign(uYcore(can_var,:)*Ycore(can_var,:,tr)') == -1
                Ycore(can_var,:,tr) = -Ycore(can_var,:,tr);
            end
        end
    end
    uYcore = mean(Ycore,3);
    sdYcore = std(Ycore,1,3);
    for e = 1:size(mB,2)
        figure(fig_ycore(e))        
        bar(uYcore(e,:));
        hold on
        h = errorbar(uYcore(e,:),sdYcore(e,:));        
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
        title(['Canonical Loadings:' num2str(e)]);
        ylabel('Correlation');
        hold off
        drawnow
    end

    [Xcore(:,:,tr) XcoreP] = corr(mU,X(TrainMaster,:));
    if tr > 1
        uXcore = mean(Xcore(:,:,1:end-1),3);
        for can_var = 1:size(uXcore,1)
            if sign(uXcore(can_var,:)*Xcore(can_var,:,tr)') == -1
                Xcore(can_var,:,tr) = -Xcore(can_var,:,tr);
            end
        end
    end
    uXcore = mean(Xcore,3);
    sdXcore = std(Xcore,1,3);
    for e = 1:size(mA,2)
        figure(fig_xcore(e));
        bar(uXcore(e,:));
        hold on
        h = errorbar(uXcore(e,:),sdXcore(e,:));
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',[LAB.spec';LAB.shape(1:3)'],'XTick',1:numel([LAB.spec';LAB.shape(1:3)']))
        title(['Canonical Loadings:' num2str(e)]);
        ylabel('Correlation');
        hold off
    end


    for e = 1:size(mA,2)
        figure(over_all(e))
        plot(mU(:,e),mV(:,e),'.')
        hold on
        title(['Overall:' num2str(e)]);
    end

end

csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/compositionANDspec_gravi_average_coreX.csv',mean(Xcore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/compositionANDspec_gravi_std_coreX.csv',std(Xcore,0,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/compositionANDspec_gravi_average_coreY.csv',mean(Ycore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/compositionANDspec_gravi_std_coreY.csv',std(Ycore,0,3));
%% F) HOLD OUT for [composition] -> gravi
T = f;
close all
X = [T.predictions];
Y = [T.LAT];
for e = 1:min(size(X,2),size(Y,2))
    fig_ycore(e) = figure;
    fig_xcore(e) = figure;
    over_all(e) = figure;
end


perDraw = .5;
Ycore = [];
Xcore = [];
for tr = 1:100
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(X,1),perDraw);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(X(TrainMaster,:),Y(TrainMaster,:));
    [Ycore(:,:,tr) YcoreP] = corr(mV,Y(TrainMaster,:));
    if tr > 1
        uYcore = mean(Ycore(:,:,1:end-1),3);
        for can_var = 1:size(uYcore,1)
            if sign(uYcore(can_var,:)*Ycore(can_var,:,tr)') == -1
                Ycore(can_var,:,tr) = -Ycore(can_var,:,tr);
            end
        end
    end
    uYcore = mean(Ycore,3);
    sdYcore = std(Ycore,1,3);
    for e = 1:size(mB,2)
        figure(fig_ycore(e))        
        bar(uYcore(e,:));
        hold on
        h = errorbar(uYcore(e,:),sdYcore(e,:));        
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
        title(['Canonical Loadings:' num2str(e)]);
        ylabel('Correlation');
        hold off
        drawnow
    end

    [Xcore(:,:,tr) XcoreP] = corr(mU,X(TrainMaster,:));
    if tr > 1
        uXcore = mean(Xcore(:,:,1:end-1),3);
        for can_var = 1:size(uXcore,1)
            if sign(uXcore(can_var,:)*Xcore(can_var,:,tr)') == -1
                Xcore(can_var,:,tr) = -Xcore(can_var,:,tr);
            end
        end
    end
    uXcore = mean(Xcore,3);
    sdXcore = std(Xcore,1,3);
    for e = 1:size(mA,2)
        figure(fig_xcore(e));
        bar(uXcore(e,:));
        hold on
        h = errorbar(uXcore(e,:),sdXcore(e,:));
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',[LAB.spec],'XTick',1:numel([LAB.spec]))
        title(['Canonical Loadings:' num2str(e)]);
        ylabel('Correlation');
        hold off
    end


    for e = 1:size(mA,2)
        figure(over_all(e))
        plot(mU(:,e),mV(:,e),'.')
        hold on
        title(['Overall:' num2str(e)]);
    end

end

csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_average_coreX.csv',mean(Xcore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_std_coreX.csv',std(Xcore,0,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_average_coreY.csv',mean(Ycore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_std_coreY.csv',std(Ycore,0,3));
%% F) HOLD OUT for [composition] -> gravi
T = f;
close all
X = [T.shapeDataL(:,1:3)];
Y = [T.LAT];
for e = 1:min(size(X,2),size(Y,2))
    fig_ycore(e) = figure;
    fig_xcore(e) = figure;
    over_all(e) = figure;
end


perDraw = .5;
Ycore = [];
Xcore = [];
for tr = 1:100
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(X,1),perDraw);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(X(TrainMaster,:),Y(TrainMaster,:));
    [Ycore(:,:,tr) YcoreP] = corr(mV,Y(TrainMaster,:));
    if tr > 1
        uYcore = mean(Ycore(:,:,1:end-1),3);
        for can_var = 1:size(uYcore,1)
            if sign(uYcore(can_var,:)*Ycore(can_var,:,tr)') == -1
                Ycore(can_var,:,tr) = -Ycore(can_var,:,tr);
            end
        end
    end
    uYcore = mean(Ycore,3);
    sdYcore = std(Ycore,1,3);
    for e = 1:size(mB,2)
        figure(fig_ycore(e))        
        bar(uYcore(e,:));
        hold on
        h = errorbar(uYcore(e,:),sdYcore(e,:));        
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',LAB.gravi,'XTick',1:numel(LAB.gravi))
        title(['Core Shape:' num2str(e)]);
        hold off
        drawnow
    end

    [Xcore(:,:,tr) XcoreP] = corr(mU,X(TrainMaster,:));
    if tr > 1
        uXcore = mean(Xcore(:,:,1:end-1),3);
        for can_var = 1:size(uXcore,1)
            if sign(uXcore(can_var,:)*Xcore(can_var,:,tr)') == -1
                Xcore(can_var,:,tr) = -Xcore(can_var,:,tr);
            end
        end
    end
    uXcore = mean(Xcore,3);
    sdXcore = std(Xcore,1,3);
    for e = 1:size(mA,2)
        figure(fig_xcore(e));
        bar(uXcore(e,:));
        hold on
        h = errorbar(uXcore(e,:),sdXcore(e,:));
        set(h,'LineStyle','none');
        set(gca,'XTickLabel',[LAB.shape],'XTick',1:numel([LAB.shape]))
        title(['Core Spec:' num2str(e)]);
        hold off
    end


    for e = 1:size(mA,2)
        figure(over_all(e))
        plot(mU(:,e),mV(:,e),'.')
        hold on
        title(['Overall:' num2str(e)]);
    end

end

csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_average_coreX.csv',mean(Xcore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_std_coreX.csv',std(Xcore,0,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_average_coreY.csv',mean(Ycore,3));
csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/composition_gravi_std_coreY.csv',std(Ycore,0,3))
%% position
for k = 1:size(f.predictions,2)
    UQ = unique(f.genoType);
    for u = 1:numel(UQ)
        fidx = find(strcmp(UQ{u},f.genoType));
        tmp = f.predictions(fidx,k);
        tmp = tmp - mean(tmp);
        MASS(fidx)  = tmp;
    end
    K = mean(f.predictions(:,k));
    K=1;
    [r,m,b] = regression(f.xpos',MASS/K);
    plot(f.xpos,MASS/K,'.');
    title(LAB.spec{k});
    waitforbuttonpress
end
%% position and growth
for k = 1:size(f.LAT,2)
    UQ = unique(f.genoType);
    for u = 1:numel(UQ)
        fidx = find(strcmp(UQ{u},f.genoType));
        tmp = f.LAT(fidx,k);
        tmp = tmp - mean(tmp);
        MASS(fidx)  = tmp;
    end
    %K = mean(f.predictions(:,k));
    [r,m,b] = regression(f.xpos',MASS);
    plot(f.xpos,MASS,'.');
    r
    waitforbuttonpress
end



%% save spec data and shape data
csvwrite(['/mnt/scratch1/phytoM/flashProjects/maize/shape_spec.csv'],[nA.specData nA.SHAPE cell2mat(nA.kernel_id)]);
fullD = [nA.specData nA.SHAPE cell2mat(nA.kernel_id)];
save(['/mnt/scratch1/phytoM/flashProjects/maize/shape_spec.mat'],'fullD','format','-double');