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
ridx = find(any(abs(diff(f.tipAngle,1,2)) > 10*pi/180,2) | ...
            any(isnan(f.specData),2) | any(isnan(f.tipAngle),2) | ...
            any(abs(f.tipAngle(:,1)) > 45*pi/180,2) | ...
            any(abs(f.tipAngle(:,1) - f.tipAngle(:,end)) < 10*pi/180,2) | ...
            any(f.shapeData(:,4) > 3*10^5,2) | ...
            any(f.shapeData(:,5) > 400,2)  | ...
            any(abs(diff(f.length,1,2)) > 4,2)  | ...
            any(f.shapeDataL(:,end) > .5*10^9) | ...
            f.predictions(:,5) > f.predictions(:,6) | ... 
            isnan(f.xpos) ...
            );  
f.tipAngle(ridx,:) = [];
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
LAB.shapeOLD = fieldString(1005:1008);
LAB.shape = {'Height' 'Width' 'Depth' 'Volume'};
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
M.length = [f.length;Mf.length];
M.LAT = [f.LAT;Mf.LAT];
M.GR = [f.GR;Mf.GR];
M.shapeDataL = [f.shapeDataL;Mf.shapeDataL];
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



%% BEGIN ANALYSIS
%% A) no hold out for [composition] to [shape]
T = f;
close all
X = [T.predictions(:,1:end)];
Y = [T.shapeDataL(:,1:3)];

%{
X = zscore(X);
Y = zscore(Y);

[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,size(X,2));
[S Y U Ey L ERR LAM] = PCA_FIT_FULL(Y,size(Y,2));
%}
%{
for e = 1:size(Ex,2)
    figure;
    bar(Ex(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

% PCA decomposition of seed shape
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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end


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
%% A) FIGURE 2 for paper
T = f;
close all
X = [T.specData];
Y = [T.shapeDataL(:,1:3)];

Y = zscore(Y);

[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);

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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end

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
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end
%% 

L = {'Weight' 'NMR' 'mg oil' 'oil percent'	'protein per' 'mg protein' 'total density' 'material density' 'total volume' 'material volume' 'air space'};
% biplot

c1 = [2^-.5 2^-.5]*[mU(:,1) mV(:,1)]';
c2 = [2^-.5 2^-.5]*[mU(:,2) mV(:,2)]';
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

%plot(mU(:,1)*norm(mA(:,1)),mU(:,2)*norm(mA(:,2)),'.');
plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');

[C_lessP] = PCA_REPROJ(lessP,Ex,Ux);
rm = [2];
L(rm) = [];
C_lessP(rm,:) = [];


% try to double check
for e = 1:size(mU,2)
    for l = 1:size(f.predictions,2)
        CHECK(e,l)= corr(f.predictions(:,l),mU(:,e));
    end
end
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,2)
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    %e1 = mA(:,1)/norm(mA(:,1));
    %p1 = e1'*C_lessP(e,:)';
    %e2 = mA(:,2)/norm(mA(:,2));
    %p2 = e2'*C_lessP(e,:)';
    %P = [p1 p2];    
    P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
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
        CHECK(e,l)= corr(Y(:,l),mU(:,e));
    end
end

L = {'Height' 'Width' 'Depth'};
C_lessP = CHECK;
for e = 1:size(C_lessP,2)    
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    %e1 = mA(:,1)/norm(mA(:,1));
    %p1 = e1'*C_lessP(e,:)';
    %e2 = mA(:,2)/norm(mA(:,2));
    %p2 = e2'*C_lessP(e,:)';
    %P = [p1 p2];    
    P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end

%{
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = mU(fidx,1:2);
    if size(sub,1) > 20
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        [tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
    end
end
%}
%% A) FIGURE 2 for paper
T = f;
close all
X = [T.specData];
Y = [T.shapeDataL(:,1:3)];

Y = zscore(Y);

[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end


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
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end
%%
L = {'Weight' 'NMR' 'mg oil' 'oil percent'	'protein per' 'mg protein' 'total density' 'material density' 'total volume' 'material volume' 'air space'};
% biplot
c1 = [2^-.5 2^-.5]*[mU(:,1) mV(:,1)]';
c2 = [2^-.5 2^-.5]*[mU(:,2) mV(:,2)]';
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
% cv
plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');



[C_lessP] = PCA_REPROJ(lessP,Ex,Ux);
rm = [2];
L(rm) = [];
C_lessP(rm,:) = [];



textMAG = 5;
MAG = 100;
hold all;
for e = 1:size(C_lessP,1)
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    e1 = mA(:,1)/norm(mA(:,1));
    p1 = e1'*C_lessP(e,:)';
    e2 = mA(:,2)/norm(mA(:,2));
    p2 = e2'*C_lessP(e,:)';
    P = [p1 p2]/MAG;    
    %P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end
%%
MAG = 1;
C_lessP = eye(3);
L = {'Height' 'Width' 'Depth'};
for e = 1:size(C_lessP,1)  
    e1 = mB(:,1)/norm(mB(:,1));
    p1 = e1'*C_lessP(e,:)';
    e2 = mB(:,2)/norm(mB(:,2));
    p2 = e2'*C_lessP(e,:)';
    P = [p1 p2];    
    %P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1)/MAG,P(2)/MAG,'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
















%% A) no hold out for [raw_spec] to [ICA_seedshape]
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];
Y = unmix';

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





%%
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];

PRE = [];
ACC = [];
T = M;
UQ = unique(T.genoType);
for tr = 1:numel(UQ);
    
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    tC = T.LAT;
    sC = [T.predictions];    
    
    TrainMaster = ~strcmp(T.genoType,UQ{tr});
    TestMaster = strcmp(T.genoType,UQ{tr});
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);    
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
    
    [gmA,gmB,gmr,gmU,gmV,gmstats] = canoncorr(sC,tC);
    gCORE = corr(gmU,T.predictions);
    
    
    PRE = [PRE;Test_sC*mA(:,:,tr)];
    ACC = [ACC;Test_tC*mB(:,:,tr)];
    
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];

end


for e = 1:size(CORE,3)
    for cc = 1:size(CORE,1)
        if sign(squeeze(CORE(cc,:,e))*gCORE(cc,:)') == -1
            CORE(cc,:,e) = -CORE(cc,:,e);
        end
    end
end
K = abs(bsxfun(@minus,CORE,gCORE));
K = squeeze(K);
K = squeeze(sum(K,2));
sCORE = sum(abs(gCORE),2);
K = bsxfun(@times,K,sCORE.^-1);
[sK sidx] = sort(K(1,:));
sUQ = UQ(sidx);
figure
bar(sK)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))
sortCORE = CORE(:,:,sidx);



% rank by 
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    mv(u,:) = sum(gmU(fidx,:));
end
[sMV sidx] = sort(mv(:,1));
sUQ = UQ(sidx);
figure
bar(sMV)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))

figure;hold all
for e = 1:size(mv,1)
    plot(abs(mv(e,1))/max(abs(mv(:,1))),abs(K(1,e))/max(abs(K(1,:))),'.');
end
legend(UQ)


%%
close all
plot(M.predictions(:,end-2),M.shapeDataL(:,end),'.')




%% try yet again
TEN = zeros(size(tC,2),size(sC,2),size(tC,1));
for e = 1:size(M.tipAngle,1)
    TEN(:,:,e) = tC(e,:)'*sC(e,:)*(sC(e,:)*sC(e,:)')^-1;
    e
end
%%
TENs = reshape(TEN,[size(TEN,1)*size(TEN,2) size(TEN,3)]);
[mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(TENs',3);
plot3(mC(:,1),mC(:,2),mC(:,3),'.')
%% 
sTENs = reshape(mS',size(TEN));
for e = 1:size(sC,1)
    tmpA = sTENs(:,:,e)*sC(e,:)';
    ang = PCA_BKPROJ(tmpA',tE,tU);
    plot(ang,'r');
    hold on;
    plot(M.tipAngle(e,:),'b');
    hold off
    drawnow
    pause(.5);
end

    










%% MASTER COMPLEX
T = M;
T = f;
T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st grouping - genotyp
T.groupingFactors(1).toGroup = 1;
T.groupingFactors(1).Groups = T.genoType;
T.groupingFactors(1).toDisplay = 1;
T.groupingFactors(1).Name = 'genoType';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd grouping - 50% whole
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
T.groupingFactors(2).toGroup = 0;
T.groupingFactors(2).Groups = holdOutGroups;
T.groupingFactors(2).toDisplay = 0;
T.groupingFactors(2).Name = 'random hold out across genotype';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3rd grouping - X percent % per genoType
UQ = unique(T.genoType);
perDrawTest = .75;
Train = zeros(size(T.specData,1),1);
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [sTrain, sTest] = crossvalind('HoldOut',size(fidx,1),perDrawTest);
    subI = find(sTrain);
    Train(fidx(subI)) = 1;
end
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
T.groupingFactors(3).toGroup = 0;
T.groupingFactors(3).Groups = holdOutGroups;
T.groupingFactors(3).toDisplay = 0;
T.groupingFactors(3).Name = 'random structured hold out accross genotype';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4th grouping kmeans
kidx = kmeans(T.specData,10);
for e = 1:numel(kidx)
    holdOutGroups{e} = num2str(kidx(e));
end
T.groupingFactors(4).toGroup = 0;
T.groupingFactors(4).Groups = holdOutGroups;
T.groupingFactors(4).toDisplay = 0;
T.groupingFactors(4).Name = 'kmeans display';


% configure pca for domain and codomain
type.PCAX.perform = 0;
type.PCAX.dim = 15;
type.PCAY.perform = 1;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
%type.learnMethod.which = 'cca';
type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
%type.learnMethod.which = 'lookup';
type.learnMethod.numComp = 15;
display.pca.perform.level = 0;
display.pca.perform.across = 1;
display.grouping = [];
display.grouping(1).grp = T.genoType;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO(T,holdOutGroups,type,display);


%% MASTER SIMPLE
T = f;
%T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));

perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
holdOutGroups = [T.genoType];


% configure pca for domain and codomain
type.PCAX.perform = 0;
type.PCAX.dim = 50;
type.PCAY.perform = 0;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
display.indivdual.perform = 0;
%type.learnMethod.which = 'cca';
%type.learnMethod.which = 'kcca';
type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
%type.learnMethod.which = 'lookup';
type.learnMethod.numComp = 15;
display.pca.perform = 1;
display.grouping = [];
display.grouping(1).grp = T.genoType;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO_OLD(T,holdOutGroups,type,display);

%% MASTER with loop over factors
cnt = 1;
RET = {};
for l = 2:1:30

    T = f;
    T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
    perDrawTest = .5;
    [Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
    for e = 1:numel(Train)
        holdOutGroups{e} = num2str(Train(e));
    end
    holdOutGroups = [T.genoType];

    % configure pca for domain and codomain
    type.PCAX.perform = 1;
    type.PCAX.dim = l;
    type.PCAY.perform = 1;
    type.PCAY.dim = 3;
    type.learnMethod.means = 0;
    type.learnMethod.which = 'cca';
    %type.learnMethod.which = 'pls';
    %type.learnMethod.which = 'net';
    %type.learnMethod.which = 'kMani';
    type.learnMethod.which = 'lookup';
    type.learnMethod.numComp = l;
    display.pca.perform = 0;
    display.indivdual.perform = 0;
    display.pca.domain = 'Spec Data';
    display.pca.codomain = 'Tip Angle';
    RET{cnt} = GO(T,holdOutGroups,type,display);
    cnt = cnt + 1;
    l
    
    
    RHO = [];
    for e = 1:numel(RET)    
        [rho pval] = corr(RET{e}.compPreGrpMeanY,RET{e}.compAccGrpMeanY);
        RHO(e,:) = diag(rho);
    end
    plot(RHO);
    drawnow
end
%% MASTER with shape
T = f;
T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
holdOutGroups = [T.genoType];
% configure pca for domain and codomain
type.PCAX.perform = 1;
type.PCAX.dim = 100;
type.PCAZ.perform = 0;
type.PCAZ.dim = 100;
type.PCAY.perform = 1;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
type.learnMethod.which = 'cca';
%type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
type.learnMethod.numComp = 60;
display.pca.perform = 1;
display.indivdual.perform = 1;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO_SUM(T,holdOutGroups,type,display);
%% simple
X = M.specData;
X = bsxfun(@minus,X,mean(X,1));
X = bsxfun(@times,X,std(X,1,1));
Y = M.tipAngle;
Y = bsxfun(@minus,Y,mean(Y,1));
Y = bsxfun(@times,Y,std(Y,1,1));
sXY = (X'*Y);
for i = 1:size(X,2)
    for j = 1:size(Y,2)
        [RHO(i,j) PVAL(i,j)] = corr(X(:,i),Y(:,j));
    end
    i
end
%% explore again
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.tipAngle,3);
for e = 1:size(xS,1)
    plot(xS(e,:),'r')
    hold on
    plot(M.tipAngle(e,:),'b')
    hold off
    drawnow
    pause(.3)
end
%%
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.specData,3);
plot3(xC(:,1),xC(:,2),xC(:,3),'.','MarkerSize',3);
hold on
CL = {'r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','g','c','r','b','k','g','c'};
for u = 1:numel(UQ)
    idx = find(strcmp(M.genoType,UQ{u}));
    if numel(idx) > 3
        %plot3(xC(idx,1),xC(idx,2),xC(idx,3),[CL{u} 'o'],'MarkerSize',5);
        [xxS xxC xxU xxE xxL xxERR xxLAM] = PCA_FIT_FULL(xC(idx,:),3);
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,1),xxE(2,1),xxE(3,1),.1,'b');
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,2),xxE(2,2),xxE(3,2),.1,'g');
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,3),xxE(2,3),xxE(3,3),.1,'r');
    end
end
%% tray var model
for u = 1:numel(UQ
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.specData,3);






%% 1) spec --> shape -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-shape/';
mkdir(oPath);
Y = f.shapeData;
X = f.specData;
dimX = 10;
dimY = 9;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);


lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Top\_Area';
L{2} = 'Top\_Major';
L{3} = 'Top\_Minor';
L{4} = 'Front\_Area';
L{5} = 'Front\_Major';
L{6} = 'Front\_Minor';
L{7} = 'Side\_Area';
L{8} = 'Side\_Major';
L{9} = 'Side\_Minor';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end
%{
%% 2) spec --> tip --  create tip angle from spec data with hold out but group over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));
%X = gradient(X);
%X = gradient(X);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));
nG = 1;
dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
X = bsxfun(@minus,X,mean(X,1));
MV = figure;
timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(f.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

[wx,wy] = myCCA(xC(Train==1,:),yC(Train==1,:),3);
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),15);

% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([20],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}
for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};



        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];

oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%}



















%% 2) NAM ONLY - PLS - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY - CCA - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - PLS - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp(G,UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
%AUG = PCA_BKPROJ(BETA',xE,xU);
%inv(B)
%figure;plot(AUG');
%tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(X,1);
    dC = mean(X,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,X,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = uSPEC;
    %C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - CCA - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end


for u = 1:numel(UQ)
    
    fidx = find(strcmp(G,UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end


%% 2) NAM ONLY PLS - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);

for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}));
    fidxTR = find(~strcmp([f.genoType],UQ{u}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY CCA - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = f.genoType;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);

for u = 1:numel(UQ)
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    [A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM CCA - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = M.genoType;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);



timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);

for u = 1:numel(UQ)
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    [A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end


%% 2) NAM ONLY - CCA - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = [f.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 5;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

% try clustering the X
kidx = kmeans(X,28);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    [A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    BETA = A*inv(B);
    
    fidx = find(strcmp([f.genoType],UQ{u}));
    fidxTR = find(~strcmp([f.genoType],UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY - PLS - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = [f.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 400;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    %BETA = A*inv(B);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(uGCx,uGCy,15);
    BETA = BETA(2:end,:);
    
    fidx = find(strcmp(G,UQ{u}));
    %fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
    %% 5) spec --> MUT tip --  predict the classic muts from the NAM parents MEANS
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
xC = PCA_REPROJ(X,xE,xU);
yC = PCA_REPROJ(Y,yE,yU);
UQ = unique(Mf.genoType);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;



for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);


%% 2) CCM & NAM - CCA - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 5;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    [A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    BETA = A*inv(B);
    
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
figure;plot(AUG');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - PLS - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 400;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    %BETA = A*inv(B);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(uGCx,uGCy,15);
    BETA = BETA(2:end,:);
    
    fidx = find(strcmp(G,UQ{u}));
    %fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end




%% 2) PLS - spec --> final_tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end
    %% 2) PLS - spec --> MUT final_tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);




Pre = [];
Acc = [];
uCPre = [];
uCAcc = [];


UQ = unique(Mf.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
%perDraw = .5;
%[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(Y(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


%% 2) PLS - spec --> avg_swingRate --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end
    %% 2) PLS - spec --> MUT avg_swingRate --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);




Pre = [];
Acc = [];
uCPre = [];
uCAcc = [];


UQ = unique(Mf.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
%perDraw = .5;
%[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(Y(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


%%
%%%%%%
%%


%%  2.25) spec --> final angle -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-finalTipAngle/';
mkdir(oPath);
Y = f.tipAngle(:,end);
X = f.specData;
dimX = 10;
dimY = 1;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Final Tip Angle';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end
%%  2.35) spec --> average swing rate -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-averageSwingRate/';
mkdir(oPath);
Y = mean(diff(f.tipAngle,1,2),2);
X = f.specData;
dimX = 10;
dimY = 1;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Average Swing Rate';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end



%%      2.5) spec --> tip angle robust looped sample for predictions
Y = f.tipAngle;
Y = bsxfun(@minus,Y(:,1),Y);
X = f.specData;
dimX = 100;
dimY = 3;
scale = 180/pi;
nG = 1;
cor = figure;
pv = figure;
R = [];
pV = [];
for L = 1:50
    perDrawTest = .5;
    [Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
    [xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
    [yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
    
    
    %net = feedforwardnet(20);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    lM = lManifold();
    lM.setmodelCompX(dimX);
    lM.setmodelCompY(dimY);
    lM.addXY(xC(Train==1,:),yC(Train==1,:));
    lM.setGroupN(nG);
    lM.learn();
    
    %subY = net(xC(Test==1,:)')';
    subY = lM.predict(xC(Test==1,:));
    
    
    subY = PCA_BKPROJ(subY,yE,yU);
    subAY = Y(Test==1,:);
    [RHO,PVAL] = corr(subY,subAY);
    figure(cor);
    hold on
    plot(diag(RHO));
    R = [R diag(RHO)];
    figure(pv);
    hold on
    plot(diag(PVAL));
    pV = [pV diag(PVAL)];
end
figure(cor);
plot(mean(R,2),'k','LineWidth',5);
figure(pv);
plot(mean(pV,2),'k','LineWidth',5);

% display the groups in doamin
figure;
LEG = {};
for e = 1:size(lM.Ux,1)
    centerX = PCA_BKPROJ(lM.Ux(e,:),xE,xU);
    plot(centerX);
    LEG{e} = num2str(e);
    hold all
end
plot(mean(X),'k');
legend(LEG)
%p1 = anova1(lM.rawX(:,2),num2str(lM.kidx));
figure;
LEG = {};
for e = 1:size(lM.Uy,1)
    kidx = (lM.kidx == e);
    centerY = PCA_BKPROJ(lM.rawY(kidx,:),yE,yU);
        LEG{e} = num2str(e);
    ugY = mean(centerY,1);
    sgY = std(centerY,1,1)*size(centerY,1).^-.5;
    errorbar(ugY,sgY);
    hold all
end
legend(LEG);
%p1 = anova1(lM.rawY(:,1),num2str(lM.kidx));
%p2 = anova1(lM.rawY(:,2),num2str(lM.kidx));
%p3 = anova1(lM.rawY(:,3),num2str(lM.kidx));
plot(mean(Y),'k');
%% 3) shape --> tip --  create tip angle from spec data with hold out but group over genotype

close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.shapeData;
Y = bsxfun(@minus,Y(:,1),Y);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 9;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = []
Acc = [];
UQ = unique(f.genoType);
scale = -180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% neural network
net = feedforwardnet(3);
net.trainParam.max_fail = 20;
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        k = lM.predict(xC(fidx,:));
        k = net(xC(fidx,:)')';
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];




        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
    end

end





for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['SHAPE --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 3) shape X spec --> tip --  create tip angle from spec data with hold out but group over genotype

close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-tip/';
mkdir(outPath);
disp = 1;


Y = f.tipAngle;
X = f.shapeData;
Z = f.specData;
Y = bsxfun(@minus,Y(:,1),Y);


% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

dimX = 8;
dimY = 3;
dimZ = 8;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[zS zC zU zE zL zERR zLAM] = PCA_FIT_FULL(Z,dimZ);
X = [xC zC];
dimX = dimX+dimZ;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
xC = X;


net = network();
net.numInputs = size(xC,2);
net.numLayers = 3;


% neural network
net = feedforwardnet(9);
net = fitnet(9);
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');


timer = .1;
YP = [];
XP = [];
Pre = []
Acc = [];
UQ = unique(f.genoType);
scale = -180/pi;


% clustered learning manifold
nG = 1;
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();




for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        k = lM.predict(xC(fidx,:));
        k = net(xC(fidx,:)')'; 
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];




        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
    end

end





for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['SHAPE --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 4) shape --> tip -- nonlinear kmeans manifold with hold out on genotype
Y = f.tipAngle;
X = f.shapeData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 10;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
PX = [];
PY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end
%% 5) spec --> MUT tip --  predict the classic muts from the NAM parents
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;
%p = randperm(size(X,1));
%Y = Y(p,:);
Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
myC = PCA_REPROJ(Mf.tipAngle,yE,yU);
mxC = PCA_REPROJ(Mf.specData,xE,xU);

timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% neural network
%net = feedforwardnet([15]);
net = fitnet([20]);
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        k = lM.predict(mxC(fidx,:));
        
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;myC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 5) PLS spec --> MUT tip --  predict the classic muts from the NAM parents
%{
tmp = f;
f = Mf;
Mf = tmp;
%}
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;


Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 100;
dimY = 3;
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
myC = PCA_REPROJ(Mf.tipAngle,yE,yU);

% simple pls regress
X = bsxfun(@minus,X,mean(X,1));
mxC = bsxfun(@minus,Mf.specData,mean(f.specData,1));
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),5);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);


timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;


%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %k = lM.predict(mxC(fidx,:));
        k = mxC(fidx,:)*BETA;
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %if strcmp(UQ{u}(end),'t')
            % gather predictions
            Pre = [Pre;k];
            Acc = [Acc;myC(fidx,:)];
        %end
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);
        
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
              %{  
                figure(MV);
                errorbar(u2*scale,s2*scale,'b');
                hold on
                errorbar(u1*scale,s1*scale,'m');
            %}

                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                hold off
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                %if mod(u,2) == 0
                %    hold off
                %end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 5) PLS spec --> MUT tip --  predict the classic muts from the NAM parents
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;
p = randperm(size(X,1));
%Y = Y(p,:);
Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 100;
dimY = 3;
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

% simple pls regress
X = bsxfun(@minus,X,mean(X,1));
mxC = bsxfun(@minus,Mf.specData,mean(f.specData,1));
UX = [];
UY = [];
UQ = unique(f.genoType);
for u = 1:numel(UQ) 
    fidx = find(strcmp([f.genoType],UQ{u}));
    UX = [UX;mean(X(fidx,:),1)];
    UY = [UY;mean(yC(fidx,:),1)];
end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),10);
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(UX,UY,10);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);


timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;


%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %k = lM.predict(mxC(fidx,:));
        k = mxC(fidx,:)*BETA;
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        if strcmp(UQ{u}(end),'t')
            % gather predictions
            Pre = [Pre;k];
            Acc = [Acc;myC(fidx,:)];
        end
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);





%% try nonlinear kmeans manifold with hold out on genotype - predict means and ind from means
Y = f.tipAngle;
X = f.specData;
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;

for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(subX);

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(3);
        
        lMStore{u} = lM;
    end
end

[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end


%% 2) PLS spec --> tip --  create tip angle from spec data with hold out but group over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));
%X = gradient(X);
%X = gradient(X);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));
nG = 1;
dimX = 100;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
X = bsxfun(@minus,X,mean(X,1));
MV = figure;
timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(f.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);


% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),15);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);

% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([20],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}
for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA;
        %k = net(xC(fidx,:)')';
        
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};



        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
%AUG = PCA_BKPROJ(BETA',xE,xU);
AUG = BETA';
figure;plot(AUG');
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    %delta = delta'*[0];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    %uSPEC = bsxfun(@plus,uSPEC,dC);
    %C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = uSPEC;
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
    axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end



%% pie graphs
for g = 1:nG
    sidx = find(Train==1);
    fidx = find(lM.kidx == g);
    gs = f.genoType(sidx(fidx));
    for u = 1:numel(UQ)
        tG(u) = sum(strcmp(gs,UQ{u}));
    end
    tG = tG / sum(tG);
    figure;
    pie(tG);
        title(['Group' num2str(g)]);
    legend(UQ);
end
%% try nonlinear kmeans manifold with hold out on genotype
Y = f.tipAngle;
X = f.specData;
%Y = bsxfun(@minus,Y,Y(:,1));
%Y = diff(Y,1,2);
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
UQ = unique(f.genoType);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        lM.addXY(xC(uidx,:),yC(uidx,:));
        lM.setGroupN(nG);
        lM.learn();
    
        %net = fitnet(20);
        %net = train(net,xC(uidx,:)',yC(uidx,:)');
        
        subX = xC(Test,:);
        subAY = Y(Test,:);
        
        subY = lM.predict(subX);
        %subY = net(subX')';
        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        [RHO,PVAL] = corr(subY,subAY);
        figure;plotyy(1:size(RHO,1),diag(RHO),1:size(RHO,1),diag(PVAL));
        title(UQ{u});
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(.1);
    end
end
[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX',uPY');
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    for u = 1:numel(UQ)
        plot(mean(cPX(pUQ==u,e)),mean(cPY(pUQ==u,e)),'r*');
    end
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title([num2str(RHO) '--' num2str(PVAL)]);
end



[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% BELOW predicitons are with genotype means for building model with hold out
%% spec --> tip -- nonlinear kmeans manifold with hold out on genotype spec
Y = f.tipAngle;
X = f.specData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end

% p
for g = 1:size(lM.pVec,1)
    pV = PCA_BKPROJ(squeeze(lM.pVec(g,:,:))',xE,xU);
    figure;
    plot(pV')
    title([num2str(g)]);
end
% look at 
for g = 1:numel(lM.corVecA)
    pV = PCA_BKPROJ(squeeze(lM.corVecA{g})',xE,xU);
    figure;
    plot(pV')
    title([num2str(g)]);
end
%{
hold on
plot(pC(1,:),'k');
pC(1,:)*lM.corVecA{g}'
%}
% 
for g = 1:numel(lM.corVecB)
    pV = PCA_BKPROJ(squeeze(lM.corVecB{g})',yE,yU);
    figure;
    plot(pV')
    title([num2str(g) '-- with means added' ]);
end
%
for g = 1:numel(lM.corVecB)
    pV = PCA_BKPROJ(squeeze(lM.corVecB{g})',yE,yU);
    pV = bsxfun(@minus,pV,yU);
    figure;
    plot(pV')
    title([num2str(g) '-- with means subtracted' ]);
end
%% spec x shape -- > tip -- nonlinear kmeans manifold with hold out on genotype 
Y = f.tipAngle;
X = f.shapeData;
Z = f.specData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 3;
dimY = 3;
dimZ = 8;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[zS zC zU zE zL zERR zLAM] = PCA_FIT_FULL(Z,dimZ);
X = [xC zC];
dimX = dimX+dimZ;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end
% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end
%% try nonlinear kmeans manifold with hold out on genotype - predict means and ind from means
Y = f.tipAngle;
X = f.specData;
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;

for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(subX);

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(3);
        
        lMStore{u} = lM;
    end
end

[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end


%% look at indivduals tip angles per genotype
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}));
    figure;
    plot(f.tipAngle(uidx,:)','g--');
    hold on
    plot(mean(f.tipAngle(uidx,:)),'k','LineWidth',5)
end
%% look at oil,protien,density,volume
X = f.tipAngle;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
for p = 1:size(f.predictions,2)
    for e = 1:size(xC,2)
        figure
        plot(f.predictions(:,p),xC(:,e),'.');
        title([num2str(p) '--' num2str(e)]);
    end
end
%% create prediction plot(s) of values for tip angle vs oil,protien, density, volume
% build holdout index
perDraw = .01;
Y = f.tipAngle;
X = f.predictions;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
plotSet = Train;
dimX = 4;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[tS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);    


%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
% predict with test
k = xC(plotSet==1,:)*A*inv(B);
for i = 1:size(yC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% kmeans on oil and plot tip angle
NG = 5;
oil_idx = kmeans(f.predictions(:,2),NG);
clear LEG
for g = 1:NG
    uidx = find(oil_idx==g);
    uT = mean(f.tipAngle(uidx,:),1);
    sT = std(f.tipAngle(uidx,:),1,1)*numel(uidx)^-.5;
    errorbar(uT,sT);
    hold all
    uP(g) = mean(f.predictions(uidx,2));
    LEG{g} = num2str(uP(g));
end
legend(LEG)




%% normalize
for e = 1:size(pC,1)
    npC(e,:) = pC(e,:)/norm(pC(e,:));
end
npC*npC'



















%% 0) load shape data
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
S = schemas(conn);
T = tables(conn);
P = columns(conn);
sid = 'kernel_3d';
for e = 1:size(P,1)
    if strcmp(P{e},sid)
        cNames = P{e,2};
        sidNUM = e;
    end
end
%% 0.1) construct query for kernel 3D
q = ['select * from kernel_3d '...    
    'join kernels on kernel_3d.kernel_id = kernels.id ' ...
    'join kernel_plates on kernels.plate_id = kernel_plates.id '...
    'join population_lines on kernel_plates.population_line_id = population_lines.id'];
results = fetch(conn,q);
%% 0.2) put kernel features into hashmap
import java.util.HashMap;
kernelVec = HashMap();
ROWnum = 6;
COLnum = 8;
ROWvec = {'A','B','C','D','E','F'};
COLvec = {'1','2','3','4','5','6','7','8'};
for e = 1:size(results,1)
    
    try
        NUMpos = results{e,45}-1;
        
        
        if NUMpos > ROWnum*COLnum
            NUMpos = NUMpos - ROWnum*COLnum;
        end
        plateN = results{e,50};
       
        if strcmp(plateN(end),'A') | strcmp(plateN(end),'B')  
            plateN = [plateN(1:end-1) '-' plateN(end)];    
        end
        
        
        
        colN = mod(NUMpos,COLnum) + 1;
        colV = ['_' COLvec{colN}];

        rowN = floor(NUMpos/COLnum) + 1;
        rowV = ['_' ROWvec{rowN}];
        
        
    catch ME

        
    end
    
    
    wellN = [rowV colV];

    key1 = [plateN '*' wellN];
    if strfind(plateN,'07S-MO001')
        plateN
    end
    %fprintf([key1 '-->' wellN '-->' num2str(NUMpos) '\n']);
    kernelVec.put(key1,results(e,3:42));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 1) load csv masterlist
[masterList, result]= readtext('/mnt/spaldingdata/florida/files/nir_master_list.csv');
%% 2) load oil strach protein prdict vectors
[pC] = readtext('/mnt/spaldingdata/nate/COEFF.csv');
pC = cell2mat(pC);
pC(:,end) = [];
for i = 1:size(pC,1)
    npC(i,:) = pC(i,:)/norm(pC(i,:));
end
%% 3) make platename -->genotype key value store from masterlist
import java.util.HashMap;
K1 = 7;
V1 = 11;
P2G = HashMap();
for e = 1:size(masterList,1)
    key = masterList{e,K1};
    value = masterList{e,V1};    
    if ~isempty(value)
        P2G.put(key,value);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 4) dig for tip angle data
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction/loganSPOOL/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% 6) open tip angles
cnt = 1;
DT = [];
WNT = {};
PLT = {};
for i = 1:numel(FileList)
    try
        [p n ext] = fileparts(FileList{i});    
        fidx = strfind(n,'_');

        plateName = n(fidx(1)+1:fidx(2)-1);
        % factor
        wellNames = n(fidx(2)+1:end);
        D = csvread(FileList{i});
        wellNames = ['_' wellNames '_'];
        fidx = strfind(wellNames,'_');
    
        if numel(fidx)-1 == size(D',2)
            DT = [DT D'];
            for i = 1:numel(fidx)-1
                WNT{cnt} = wellNames(fidx(i)+1:fidx(i+1)-1);
                PLT{cnt} = plateName;
                cnt = cnt + 1;
            end
        end
    catch ME
        ME
      
    end
end
%% 7) clean tip angles
% clean via derivative
ridx = find(any(abs(diff(DT,1,1)) > 20/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
% clean via total bend
ridx = find(any(abs(DT(1,:) - DT(end,:)) < 20/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 8) run import init from mongoDB
import phytoG.locked.Bpersist.Bos.implementations.*
import com.mongodb.*;
import java.util.Map.*;
import java.util.HashMap;
import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
import phytoG.locked.BdataObjects.BbioObjects.maize.spectraData;
oStore = OStore_mdb();
oStore.accessResource();
oStore.setCollection('specData');
%% 9) import spec and pull kernel shape data from hashmap
f.specData = [];
f.tipAngle = [];
f.kernelData = [];
f.weightData = [];
f.genoType = {};
f.wellNumber = {};
f.plateName = {};
for e = 1:numel(WNT)
    qMap = HashMap();
    pn = PLT{e};
    qMap.put('_k_m._pnode._k_m._plateName',pn);
    wn = ['_' WNT{e}(1) '_' WNT{e}(2)];
    qMap.put('_k_m._pnode._k_m._wellName',wn);
    cursor = oStore.search(qMap);
    itr = cursor.iterator();
    if itr.size()
        n = itr.next();
        n = spectraData(n);
        
        spec = n.getSpectrum();
        
        
        
        %%%%%%%
        % pull spec data from array
        specData = [];
        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        
        %%%%%%%
        % get kernel shape data
        sKey = [pn '*' wn];
        kVec = kernelVec.get(sKey);
        kernelData = [];
        if ~isempty(kVec)
            for i = 1:kVec.size()
                kernelData(i,1) = kVec(i);
            end
        end
        
        
        next_specData = specData;
        next_kernelData = kernelData;
        next_weight = n.getWeight();
        next_tipAngle =  DT(:,e);
        next_genoType = P2G.get(pn(1:end-2));
        next_wellName = wn;
        next_plateName = pn;
        
        %%%%% check for valid entries
        if ~isempty(next_specData) && ...
           ~isempty(next_kernelData) && ...
           ~isempty(next_tipAngle) && ...
           ~isempty(next_genoType) && ...
           ~isempty(next_wellName) && ...
           ~isempty(next_plateName) && ...
           ~isempty(next_weight)
           
           
            f.specData = [f.specData next_specData];
            f.tipAngle = [f.tipAngle next_tipAngle];
            f.weightData = [f.weightData str2num(char(next_weight))];
            f.kernelData = [f.kernelData next_kernelData];
            f.genoType{end+1} = next_genoType;
            f.wellNumber{end+1} = next_wellName;
            f.plateName{end+1} = next_plateName;
            fprintf(['plate: ' pn ' well: ' wn ' found! \n']);
        else
            fprintf(['plate: ' pn ' well: ' wn ' incomplete Data found! \n']);   
        end
        
        
        
    else
        fprintf(['plate: ' pn ' well: ' wn ' NOT found! \n']);
    end
end
f.specData(1:3,:) = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 10) get subset of kernel feature data - area and lengths
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
%% 10.5) CHOICE - look at histogram(s)
for i = 1:size(f.kernelSub,1)
    figure
    hist(f.kernelSub(i,:))
end
%% 11) clean data
% clean on area
ridx = find(f.kernelSub(1,:) > 4*10^5);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);

ridx = find(f.kernelSub(2,:) > 500);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
%% 12) get subsubset of kernel data - the lengths
sIDX = [11 30 31] - 2;
f.kernelLengths = f.kernelData(sIDX,:);
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% JUNK
disp = 1;
phytoK = {};
ridx = [];
for e = 1:numel(f.plateName)
    try
        topName = translateToKernelImage(f.plateName{e},f.wellNumber{e});
        kernelTopExtract2(topName,disp);
    end
end
%% XX-13) look at image via lookups
disp = 1;
phytoK = {};
ridx = [];
for e = 1:numel(f.plateName)
    try
        topName = translateToKernelImage(f.plateName{e},f.wellNumber{e});
        [height width area contour] = kernelTopExtract2(topName,disp);
        phytoShapes = [height width area];
        k = [f.plateName{e} '*' f.wellNumber{e}];
        phytoK(e).key = k;
        phytoK(e).value = phytoShapes;
        phytoK(e).curve = contour;
        fprintf([k '\n']);
        ridx(e) = 0;
    catch ME
        phytoK{e} = [];
        ridx(e) = 1;
        fprintf(['error@' topName '\n']);
    end
end
%% JUNK
ridx = find(ridx);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
sIDX = [11 30 31] - 2;
f.kernelLengths = f.kernelData(sIDX,:);
%% XX-13) do inserts for phytoKernel
f.phytoKernel = [];
for e = 1:numel(f.plateName)
     k = [f.plateName{e} '*' f.wellNumber{e}];
     fidx = find(strcmp({phytoK.key},k));
     f.phytoKernel(e,:) = phytoK(e).value;
end
f.phytoKernel = f.phytoKernel';
%% BREAK
%% look at lengths by lengths
close all
for i = 1:size(f.kernelLengths,1)
    for j = 1:size(f.kernelLengths,1)
        figure;
        plot(f.kernelLengths(i,:),f.kernelLengths(j,:),'.')
        title([num2str(i) '--' num2str(j)]);
    end
end
%% 11.2) create weights from lengths with hold out
close all
disp = 1;
X = f.kernelLengths';
Y = f.weightData';
perDraw = .7;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);

[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,3);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,1);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
[mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),1);
% predict with test
k = sC(Test==1,:)*A*inv(B);
mk = sC(Test==1,:)*mA*inv(mB);


k = PCA_BKPROJ(k,tE,tU);
mk = PCA_BKPROJ(mk,tE,tU);


plot(k,Y(Test==1),'.')
figure;
plot(mk,Y(Test==1),'.');

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% MUST) flip tip angle
f.tipAngle = -f.tipAngle;





%% 10) hold out
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',6);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end
%%%%%%%%
%% BREAK
%% BREAK
%% BREAK
%% lets look at volume predictions
close all
offset = 0;
plot3(f.kernelSub(offset+1,:),f.kernelSub(offset+2,:),f.kernelSub(offset+3,:),'.');
figure;
plot(f.kernelSub(offset+1,:),f.kernelSub(offset+2,:).*f.kernelSub(offset+3,:),'.');
figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.kernelSub(offset+1:offset+3,:)',3);
plot3(sC(:,1),sC(:,2),sC(:,3),'.');
%% 10.5) look at corr between weights and area,lengths and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelSub';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};

for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
end
legend(UQ);
%% 10.5) look at corr between weights and length and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelLengths';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'along length of cob - front minor','kernel height - top major','along circ of cob - top minor'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Lengths');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 125;
%% 10.5) look at corr between weights and phytoLengths and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.phytoKernel';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'along length of cob - front minor','kernel height - top major','along circ of cob - top minor'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Lengths');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 125;
%% 10.5) look at corr between weights and area and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelData(27,:)';  % only look at top area

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'top area'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Area');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(1,:) < 125;
%% 10.6) look at volumes and density vs weight(s)
% volume(S)
Y = f.weightData';
X = [f.kernelData(8,:).*f.kernelData(9,:).*.5.*(f.kernelData(18,:)+f.kernelData(28,:));...
     f.kernelData(18,:).*f.kernelData(19,:).*.5.*(f.kernelData(8,:)+f.kernelData(29,:));...
     f.kernelData(28,:).*f.kernelData(29,:).*.5.*(f.kernelData(19,:)+f.kernelData(9,:))]';
plot3(X(:,1),X(:,2),X(:,3),'r.');
title('Test each volumes vs themselves');
axis equal
%X = mean(X,2).^(1/3);
% calc density from volume and mass
X = mean(X,2);
X = [X Y.*X.^-1];

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};

for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
end
%% BREAK
%% 11) create tip angle from spec data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
for L = 1:1
    % loop over number of basis vectors
    for num = 15

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %{
        %%%%%%%%%%%%%%%
        % pls regression        
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from spec data AND shape data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = [f.specData;f.kernelSub]';
for L = 1:1
    for num = 15

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}

        %{
        %%%%%%%%%%%%%%%
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from kernel-Lengths with NO hold out
%close all
disp = 1;
Y = f.tipAngle';
X = [f.kernelLengths]';
for L = 1:1
    for num = 3

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}

        %{
        %%%%%%%%%%%%%%%
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from spec data AND shape data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = [f.specData]';
%X2 = [f.specData;f.kernelSub]';
X2 = [f.specData]';
for L = 1:1
    for num = 5

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [xS2 sC2 sU2 sE2 sL2 sERR2 sLAM2] = PCA_FIT_FULL(X2,num);
        [xS2s sC2s sU2s sE2s sL2s sERR2s sLAM2s] = PCA_FIT_FULL(f.kernelSub',6);        
        sC2 = [sC2 sC2s];
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        [A2,B2,r2,U2,V2,stats2] = canoncorr(sC2,tC);
        % predict
        k = sC*A*inv(B);
        k2 = sC2*A2*inv(B2);


        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        k2 = PCA_BKPROJ(k2,tE,tU);
        
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            % find the genotype
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u12 = mean(k2(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s12 = std(k2(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);

            dif(u) = mean(u1*u2');
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    hold on
                    errorbar(u12*-180/pi,s12*-180/pi,'g');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            %}
        end
        DIF(num) = mean(dif);
    end
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.2) create tip angle from spec data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_spec/';
disp = 1;
Y = f.tipAngle';
X = f.specData';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,15);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.2) create tip angle from shape data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_shape/';
disp = 1;
Y = f.tipAngle';
X = f.kernelLengths';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,3);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),3);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.2) create tip angle from shape data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_shapeANDspec/';
disp = 1;
Y = f.tipAngle';
X = f.kernelLengths';
Z = f.specData';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL([X Z],7);    
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

    
    
    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),3);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.3) create hold out via genotype
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
UQ = unique([f.genoType]);
num = 15;
for L = 1:3
    for u = 1:numel(UQ)

        Test = strcmp([f.genoType],UQ{u});
        Train = ~ Test;


        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train,:),tC(Train,:));
        % predict with test
        k = sC(find(Test==1),:)*A*inv(B);


        k = PCA_BKPROJ(k,tE,tU);

        u1 = mean(k);
        u2 = mean(Y(Test==1,:),1);

        s1 = std(k,1,1);
        s2 = std(Y(Test==1,:),1,1);

        
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                pause(1);
                hold off
            catch

            end
        end
    end
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.7) create prediction plot(s) of values for tip angle vs specdata
% build holdout index
perDraw = .01;
Y = f.tipAngle';
X = f.kernelSub';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

num = 6;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
%[pA,pB] = myPLS1(sC(Train==1,:),tC(Train==1,:));
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for tip angle vs shapedata
% build holdout index
perDraw = .01;
Y = f.tipAngle';
X = f.specData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for tip angle vs kernel Lengths
% build holdout index
close all
perDraw = .01;
Y = f.tipAngle';
X = f.kernelLengths';
numX = 3;
numY = 3;
plotSet = Train;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,numX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,numY);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%[A,B,U,V] = myCCA(xC(Train==1,:),yC(Train==1,:),5);
% predict with test
k = xC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = xC(plotSet==1,i);
    %yP = k(:,i);
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo');
    %axis equal
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    %plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    %axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for spec vs shapedata
% build holdout index
perDraw = .01;
Y = f.kernelSub';
X = f.specData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

numX = 6;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,numX);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for weight vs shapedata
% build holdout index
perDraw = .01;
numX = 6;
numY = 1;
X = f.kernelSub';
Y = f.weightData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,numX);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,numY);

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.8) create oil, starch, protien from tip angles via CCA
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;    
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via CCA by genotype
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);    
    figure;
    hold on
    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        realV = mean(f.ospw(i,fidx));
        preV = mean(k(fidx,i));        
        plot(realV,preV,'.');        
    end
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via pls
% build holdout index
perDraw = .3;
f.ospw = pC*f.specData;
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

% pls regression
[XL,YL,XS,YS,BETA] = plsregress(tC,sC,5);
% predict
k = [ones(size(tC,1),1) tC]*BETA;


k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;        
    %plot(f.ospw(i,Test==1),k(:,i),'.')
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.9) create oil, starch, protien from tip angles directly
% build holdout index
f.ospw = pC*f.specData;
X = f.ospw';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,4);


%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);

TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);        
    figure;
    hold on
    
    realV = (f.ospw(i,:));
    preV = (k(:,i));        
    plot(realV,preV,'.');        

    title(TITLE{i})
    axis equal
end
%% test myCCA
Y = f.tipAngle';
X = f.specData';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);
[A,B,r,U,V,stats] = canoncorr(sC,tC);
[wx,wy] = myCCA(sC,tC,3);
%%

close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
fitX = f.specData';
fitY = f.tipAngle';
cnt = 1;
for num = 3:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        
        
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(fitX,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(fitY,3);
        
        
    %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(Train==1,yi),tmpX(Train==1,:));
            k = [k tmpX*b];
        end
        %}
        


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
        % predict
        %predictY = sC*A*inv(B);
        predictY = sC*A*B';



        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}
        
        predictY = PCA_BKPROJ(predictY,tE,tU);
        
        
        trTOP = sum(predictY(Train==1,:).*fitY(Train==1,:),2);
        trBOTX = sum(predictY(Train==1,:).*predictY(Train==1,:),2);
        trBOTY = sum(fitY(Train==1,:).*fitY(Train==1,:),2);
        trCORR = trTOP.*trBOTX.^-.5.*trBOTY.^-.5;
        
        tsTOP = sum(predictY(Test==1,:).*fitY(Test==1,:),2);
        tsBOTX = sum(predictY(Test==1,:).*predictY(Test==1,:),2);
        tsBOTY = sum(fitY(Test==1,:).*fitY(Test==1,:),2);
        tsCORR = tsTOP.*tsBOTX.^-.5.*tsBOTY.^-.5;
        
        
        trainR(loop) = mean(trCORR);
        testR(loop) = mean(tsCORR);        
        
    end
    
    trU(cnt) = mean(trainR);
    trS(cnt) = std(trainR);
    
    tsU(cnt) = mean(testR);
    txS(cnt) = std(testR);
    
    cnt = cnt + 1;
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end
%% try correlation via genotype
close all
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',5);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);


UQ = unique([f.genoType]);
tipFig = figure;
specFig = figure;
ccFig = figure;
ccFigT = figure;
ccFigS = figure;
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'}
for u = 1:numel(UQ)
    -180/pi
    
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData(:,fidx)',8);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle(:,fidx)',8);
    
    
    % get the tip angle
    tmpTip = f.tipAngle(:,fidx)';
    % get the spec data
    tmpSpec = f.specData(:,fidx)';
    % get the means
    tU = mean(tmpTip,1);
    sU = mean(tmpSpec,1);
    % get the std
    tS = std(tmpTip,1,1)*size(tmpTip,1)^-.5;
    xS = std(tmpSpec,1,1)*size(tmpSpec,1)^-.5;
    
    if numel(fidx) > 1
        %[A,B,r,U,V,stats] = canoncorr(sC(fidx,:),tC(fidx,:));
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        
        
        
        
        figure(ccFigS);
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        plot(sData);
        hold all;
        
        figure(ccFigT);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(tData);
        hold all;
        
        
        LEG{u} = [UQ{u} '--' num2str(r(1)) '--' num2str(numel(fidx)) '--' num2str(stats.p(1))];
        figure(ccFig);
        scatter(U(:,1),V(:,1),CL{u});axis equal
        hold all
    else
        figure(ccFigT);
        scatter(0,0,CL{u});
        LEG{u} = [UQ{u} '---NA'];
    end
    
    
    figure(tipFig);
    errorbar(tU,tS);
    hold all
    
    figure(specFig);
    errorbar(sU,xS);
    hold all
    
    
end
figure(tipFig);
legend(UQ);
figure(specFig);
legend(UQ);
figure(ccFig);
legend(LEG);
figure(ccFigS);
legend(UQ);
figure(ccFigT);
legend(UQ);
%% try genotype cluster

f.ospw = pC*f.specData;


figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(sC(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,stats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);

scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
f.ospw = pC*f.specData;
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
%X = f.ospw';
for num = 15
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(Y(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(Y(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}%% oil starch protien

f.ospw = pC*f.specData;

figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(f.ospw(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,st
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
endats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);

scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
%% please try
f.ospw = pC*f.specData;
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
%X = f.ospw';
for num = 15
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(Y(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(Y(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}
%% size matters
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',size(f.specData',2));
newData = (diag(diag(sLAM).^-.5)*sC')';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(newData,5);
%%
ng = 1;
kidx = kmeans(f.tipAngle',ng);
kidx = kmeans(f.specData',ng);
%% init first tip angle to zero
f.tipAngle2 = bsxfun(@minus,f.tipAngle,f.tipAngle(1,:));
%f.tipAngle2 = diff(f.tipAngle,1,1);
%% 
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',8);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',8);
T = [sC tC];
ng = 1;
kidx = kmeans(T,ng);C = PCA_REPROJ(,Ey,Uy);
%% try upbulk 2 term
for i = 1:size(sC,2)
    for j = i:size(sC,2)
        sC = [sC sC(:,i).*sC(:,j)];
    end
end
sE%% explore nonlinear
plot3(sC(:,1),sC(:,2),tC(:,2),'.')
%% pls regression
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',11);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
[XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
Y0 = [ones(size(sC,1),1) sC]*BETA;
for e = 1:size(Y0,1)
    plot(Y0(e,:),'r');
    sData = PCA_BKPROJ(A(:,1)',sE,sU);
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% clustered cannon corr --> 
CL = {'r.','b.','g.','k.','m.'};
figure;
for num = 10
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
    for e = 1:ng
        [A,B,r,U,V,stats] = canoncorr(sC(kidx==e,:),tC(kidx==e,:));
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(U(:,1),V(:,1),CL{e});axis equal
        hold on;
        R(num) = r(1);
    end
    hold off;
end

axis([-4 4 -4 4])
figure;plot(sData);
figure;plot(tData);
figure;plot(-f.tipAngle);
figure;plot(f.specData);
%%
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end


%%
pth = '/mnt/spaldingdata/Bessie/forPoster/';
csvwrite([pth 'tipAngle.csv'],f.tipAngle);
%%
csvwrite([pth 'corr.csv'],[U(:,1),V(:,1)]);
%% view via eye
OP = tData;
mea = OP*f.tipAngle;
[J sidx] = sort(mea);

hold on;
for e = 1:numel(sidx)
    plot(tS','b');
    hold on
    plot(tS(sidx(e),:),'r');
    hold off
    drawnow
end
%% predict spec
A = sC'/tC';
sim = A*tC';
M = PCA_BKPROJ(sim',sE,sU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% predict tip angle
A = tC'/sC';
sim = A*sC';
M = PCA_BKPROJ(sim',tE,tU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.tipAngle(:,e),'b');
    hold off
    pause(.3)
end
    %%
    cnt = 1;
    while (itr.hasNext())
        n = itr.next();
        n = spectraData(n);
        spec = n.getSpectrum();

        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        DS(:,cnt) = specData;
        WNS{cnt} = n.getPlateName();
        PNS{cnt} = n.getWellName();
        cnt = cnt + 1
    end
%% joy again
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(f.specData,tC,30);
%[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(xsC,xtC,3);
Y = [ones(size(tC,1),1) f.specData]*BETA(1:end,:);
%Y = [ones(size(tC,1),1) xsC]*BETA(1:end,:);
plot(tC(:,1),Y(:,1),'.');
%plot(XS(:,1),YS(:,1),'.')
%%
  [A,B,r,U,V,stats] = canoncorr(f.specData,f.tipAngle);

  [A,B,r,U,V,stats] = canoncorr(specSTACK,tipSTACK);
    %[A,B,r,U,V,stats] = canoncorr(xC,yC);

figure;
plot(U(:,1),V(:,1),'.');
title([num2str(r(1)) '--' num2str(stats.p(1))]);
        % predict with test
        k = subX(Test==1,:)*A*inv(B);
%%
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData,20);
UQ = unique(f.genoType);
tipSTACK = [];
specSTACK = [];
G = {};
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType)==1);
    if numel(fidx) > 8
        [utS utC utU utE utL utERR utLAM] = PCA_FIT_FULL(tC(fidx,:),3);
        for bv = 1:size(utE,2)
            utE(:,bv) = utE(:,bv)*sign(utE(bv,bv));
            diag(utE);
        end
        [utC] = PCA_REPROJ(tC(fidx,:),utE,utU);
        utC = utC*diag(utL.^-1);
        %tipSTACK = [tipSTACK; utC*diag(utL.^-1)];
        %tipSTACK = [tipSTACK; utC];
        tipSTACK = [tipSTACK; tC(fidx,:)];
        
        
        [usS usC usU usE usL usERR usLAM] = PCA_FIT_FULL(sC(fidx,:),3);
        for bv = 1:size(usE,2)
            usE(:,bv) = usE(:,bv)*sign(usE(bv,bv));
        end
        [usC] = PCA_REPROJ(sC(fidx,:),usE,usU);
        usC = usC*diag(usL.^-1);
        %specSTACK = [specSTACK; usC*diag(usL.^-1)];
        %specSTACK = [specSTACK; usC];
        specSTACK = [specSTACK; sC(fidx,:)];
        G = [G ;f.genoType(fidx)];
    end
end
%%
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData,15);
UQ = unique(f.genoType);
tipSTACK = [];
specSTACK = [];
G = {};
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType)==1);
    if numel(fidx) > 6
        %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(sC(fidx,:),tC(fidx,:),3);
        [A,B,r,U,V,stats] = canoncorr(sC(Train,:),tC(Train,:));
        xsC(fidx,:) = U;
        xtC(fidx,:) = V;
        G = UQ{u};
    end
end
specSTACK = sxC;
tipSTACK = stC;

%%
UQ = unique(G);
DSe(1) = figure;
DSe(2) = figure;
DSe(3) = figure;

DSr(1) = figure;
DSr(2) = figure;
DSr(3) = figure;
ux = [];
uy = [];
uRx = [];
uRy = [];
[mA,mB,mr,mU,mV,mstats] = canoncorr(specSTACK,tipSTACK);

for u = 1:numel(UQ)
    Test = find(strcmp(UQ{u},G));
    Train = find(~strcmp(UQ{u},G));
    
    sUQ = setdiff(G,UQ{u});

    gX = [];
    gY = [];
    for u2 = 1:numel(sUQ)
        idx2 = find(strcmp(sUQ{u2},G));
        gX(u2,:) = mean(specSTACK(idx2,:),1);
        gY(u2,:) = mean(tipSTACK(idx2,:),1);
    end
    %Test = find(strcmp(UQ{u},f.genoType));
    %Train = find(~strcmp(UQ{u},f.genoType));
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(specSTACK(Train,:),tipSTACK(Train,:),3);
    %A = XL;
    %B = YL;
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(f.specData(Train,:),tC(Train,:),3);
    
    %plot(XS(:,1),YS(:,1),'.')
    
    [A,B,r,U,V,stats] = canoncorr(specSTACK(Train,:),tipSTACK(Train,:));
    MA(:,:,u) = A;
    MB(:,:,u) = B;
    %[A,B,r,U,V,stats] = canoncorr(gX,gY);
    u1 = mean(specSTACK(Train,:),1);
    u2 = mean(tipSTACK(Train,:),1);
    %BETA = A*inv(B);
    
    %ygit = [ones(numel(Test),1) specSTACK(Test,:)]*BETA;
    %ygit = [ones(numel(Test),1) f.specData(Test,:)]*BETA;
    %ygit = [specSTACK(Test,:)]*BETA;
    xfit = bsxfun(@minus,specSTACK(Test,:),u1)*A;
    yfit = bsxfun(@minus,tipSTACK(Test,:),u2)*B;    
    ux = [ux ;mean(xfit,1)];
    uy = [uy ;mean(yfit,1)];
    
    for e = 1:size(yfit,2)
        figure(DSe(e))
        hold on
        %plot(ygit(:,e),tC(Test,e),'b.');
        plot(xfit(:,e),yfit(:,e),'b.');        
    end
    %ygit = [specSTACK(Train,:)]*BETA(2:end,:);
    %ygit = [ones(numel(Train),1) specSTACK(Train,:)]*BETA;
    %ygit = [ones(numel(Train),1) f.specData(Train,:)]*BETA;
    %ygit = [specSTACK(Train,:)]*BETA;
    xfit = bsxfun(@minus,specSTACK(Train,:),u1)*A;
    yfit = bsxfun(@minus,tipSTACK(Train,:),u2)*B;    
    uRx = [uRx ;mean(xfit,1)];
    uRy = [uRy ;mean(yfit,1)];
    for e = 1:size(yfit,2)
        figure(DSr(e))
        hold on
        %plot(ygit(:,e),tC(Test,e),'b.');
        plot(xfit(:,e),yfit(:,e),'r.');
        plot(ux(end,e),uy(end,e),'*');
    end
end
for e = 1:size(yfit,2)
    figure(DSr(e))
    hold on
    plot(ux(:,e),uy(:,e),'k*');
    plot(uRx(:,e),uRy(:,e),'g*');
end
[RHO,PVAL] = corr(ux,uy);
%% master of joy - hold out 50% genotype of each one
T = f;
initT = T.tipAngle(:,1);
T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
UQ = unique(T.genoType);
perDraw = .5;
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,15);
%[shS shC shU shE shL shERR shLAM] = PCA_FIT_FULL(T.shapeData,3);
%sC = [sC shC];
TrainMaster = [];
TestMaster = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
    TrainMaster = [TrainMaster;Train];
    TestMaster = [TestMaster;Test];
end
TrainMaster = find(TrainMaster);
TestMaster = find(TestMaster);
u_sC = mean(sC(TrainMaster,:),1);
Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
u_tC = mean(tC(TrainMaster,:),1);
Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
[mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
Test_mU = Test_sC*mA;
Test_mV = Test_tC*mB;
BETA = mA*inv(mB);
plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end

% corr with predictions
[RHO PVAL] = corr(mU,T.predictions(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),T.tipAngle(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),T.tipAngle(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
%% joy master with mutants - predict from above
T = Mf;
initT = T.tipAngle(:,1);
T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
tC_mut = PCA_REPROJ(T.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(T.specData,sE,sU);
UQ = unique(T.genoType);
pre = figure; 
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
   
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    figure;
    errorbar(pUTA,pSTA,'r');
    hold on
    axis([0 61 0 2])
    errorbar(UTA,STA,'k');
    axis([0 61 0 2])
    hold off
    title(UQ{u});
end
%% scatter plots for mutants
initT = Mf.tipAngle(:,1);
Mf.tipAngle = bsxfun(@minus,Mf.tipAngle,initT);
tC_mut = PCA_REPROJ(Mf.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(Mf.specData,sE,sU);
UQ = unique(Mf.genoType);
preC(1) = figure; 
preC(2) = figure; 
preC(3) = figure; 
close all
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
tmpX = [];
tmpY = [];
UtmpX = [];
UtmpY = [];
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(mean(tC_mut(fidxw,e)),mean(predict(:,e)),CL{u},'MarkerSize',7);
       
       plot(linspace(-2,2,100),linspace(-2,2,100),'r');
    end
    UtmpX = [UtmpX;mean(tC_mut(fidxw,:),1)];
    UtmpY = [UtmpY;mean(predict,1)];
    
    tmpX = [tmpX;tC_mut(fidxw,:)];
    tmpY = [tmpY;predict];
    
end
%{
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(tC_mut(fidxw,e),predict(:,e),[CL{u}(1) '.']);
       plot(linspace(-10,10,100),linspace(-10,10,100),'r');
    end
end
%}
%% scatter plots for mutants
initT = Mf.tipAngle(:,1);
Mf.tipAngle = bsxfun(@minus,Mf.tipAngle,initT);
tC_mut = PCA_REPROJ(Mf.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(Mf.specData,sE,sU);
UQ = unique(Mf.genoType);
preC(1) = figure; 
preC(2) = figure; 
preC(3) = figure; 
close all
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
tmpX = [];
tmpY = [];
UtmpX = [];
UtmpY = [];
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predictSPEC = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*mA;
    predictTIP = bsxfun(@minus,tC_mut(fidxw,:),u_tC)*mB;
    
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(mean(predictSPEC(:,e)),mean(predictTIP(:,e)),CL{u},'MarkerSize',7);
       
       plot(linspace(-2,2,100),linspace(-2,2,100),'r');
    end
    
    UtmpX = [UtmpX;mean(predictSPEC,1)];
    UtmpY = [UtmpY;mean(predictTIP,1)];
    
    tmpX = [tmpX;predictSPEC];
    tmpY = [tmpY;predictTIP];
    
end
%% look at model
cv_tipAngle = PCA_BKPROJ(mB',tE,tU);
plot(cv_tipAngle(1:2,:)');
cv_specData = PCA_BKPROJ(mA',sE,sU);
non_lat_spec = predictVectors*bsxfun(@minus,cv_specData,dX)';
for e = 1:size(non_lat_spec,2)
    figure;
    bar(non_lat_spec(:,e));
end
%% try tip angle or length whole from only the predictions
T = f;%M;
%initT = T.tipAngle(:,1);
%T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
%UQ = unique(T.genoType);
perDraw = .05;
%Y = gradient(T.length);
%Y = T.length;
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
%sC = T.tot_predictions(:,[5 4 7 8 9 10 1]);
sC = T.predictions;
%sC = [T.predictions T.shapeData];
%sC = [T.predictions T.shapeData(:,[2 3 5 6 8 9])];
usC = mean(sC,1);
sC = bsxfun(@minus,sC,usC);
TrainMaster = [];
TestMaster = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
    TrainMaster = [TrainMaster;Train];
    TestMaster = [TestMaster;Test];
end
TrainMaster = find(TrainMaster);
TestMaster = find(TestMaster);
u_sC = mean(sC(TrainMaster,:),1);
Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
u_tC = mean(tC(TrainMaster,:),1);
Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
[mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
Test_mU = Test_sC*mA;
Test_mV = Test_tC*mB;
BETA = mA*inv(mB);
plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    %UTA = mean(T.tipAngle(fidxw,:),1);
    %STA = std(T.tipAngle(fidxw,:),1,1);
    
    UTA = mean(Y(fidxw,:),1);
    STA = std(Y(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end

% corr with predictions
[RHO PVAL] = corr(mU,sC(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),Y(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end

% regession along vars
for e = 1:size(mU,2)
    CS = linspace(min(mU(:,e)),max(mU(:,e)),10);
    for tm = 1:size(sC,2)
        p = polyfit(mU(:,e),sC(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(tC,2)
        p = polyfit(mV(:,e),tC(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
%% look at model factors for only predictions
close all
cv_tipAngle = PCA_BKPROJ(mB',tE,tU);
cv_tipAngle = bsxfun(@minus,cv_tipAngle,tU);
figure;
plot(cv_tipAngle');

%{
t_struct_C = (mV'*tC);
t_struct_C = PCA_BKPROJ(t_struct_C',tE,tU);
t_struct_C = bsxfun(@minus,t_struct_C,tU);
figure;
plot(t_struct_C')
%}

%{
for e = 1:size(mV,2)
    for tm = 1:size(T.tipAngle,2)
        plot(mV(:,e),T.tipAngle(:,tm),'.')
        drawnow
        pause(.1)
    end
end
%}
    
[t_struct_C pval]= corr(mV,T.tipAngle(TrainMaster,:));
figure;
plot(t_struct_C')


s_struct_C = corr(mU,T.predictions(TrainMaster,:));
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

for e = 1:size(mA,2)
    figure;
    bar(mA(:,e))
    L{e} = num2str(e);
    title(['corr_value' num2str(e)]);
end
%% lat parameters alone without hold outs
T = M;
close all
%X = [T.predictions T.shapeData];
X = [T.predictions];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,T.LAT);

[t_struct_C t_pval] = corr(mV,T.LAT);
for e = 1:size(mB,2)
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

[s_struct_C s_pval] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
end

%{
corr(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3));
figure
plot(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3),'.');
figure;
plot(mU(:,1),mV(:,1),'.');
corr(T.LAT(:,2),T.predictions(:,3));
%}

para = 1;
figure;
kidx = kmeans(mU(:,para),4);
UQ = unique(kidx);
for u = 1:numel(UQ)
    fidx = find(kidx==u);
    uT = mean(T.tipAngle(fidx,:),1);
    uS = std(T.tipAngle(fidx,:),1,1)*numel(fidx).^-.5;
    errorbar(uT,uS);
    hold all
    LEG{u} = num2str(mean(mU(fidx,para)));
end
legend(LEG)
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/lat_core_struct_tip_angle.csv',t_struct_C);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/lat_core_struct_predictions.csv',s_struct_C);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var1.csv',[mU(:,1),mV(:,1)]);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var2.csv',[mU(:,2),mV(:,2)]);
%% order NAM parents by cv1 for spec data
UQ = unique(T.genoType);
TROPICAL = [0 0 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 1 0];
%TROPICAL = [0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
TROP = [];
NON_TROP = [];
e_mU = [];
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    e_mU(u) = mean(mU(fidx));
    %{
    if TROPICAL(u)
        TROP = [TROP;f.tipAngle(fidx,:)];
    else
        NON_TROP = [NON_TROP;f.tipAngle(fidx,:)];
    end
    %}
end
[e_mU sidx] = sort(e_mU);
UQ = UQ(sidx);
close all
bar(e_mU)
%{
e_TROP = mean(TROP,1);
s_TROP = std(TROP,1,1)*size(TROP,1)^-.5;
e_NON_TROP = mean(NON_TROP,1);
s_NON_TROP = std(NON_TROP,1,1)*size(NON_TROP,1)^-.5;
errorbar(e_TROP*180/pi,s_TROP*180/pi,'r');
hold on
errorbar(e_NON_TROP*180/pi,s_NON_TROP*180/pi);
%}
%% predict shape from composition
T = f;
close all
X = [T.predictions];
Y = [T.shapeData];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
REL = 3;
[t_struct_C t_pval] = corr(mV,Y);
for e = 1:REL
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
    
end

[s_struct_C s_pval] = corr(mU,X);
for e = 1:REL
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
end
%% master of joy series - hold out 50% genotype of each one
%% 1) random holdout without respect to genotytpe
for tr = 1:100
    T = f;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    %sC = [T.predictions];
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    Test_mU = Test_sC*mA;
    Test_mV = Test_tC*mB;
    UtoV = mU\mV;
    BETA = mA*UtoV*inv(mB);

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    cCore_S_composition_m = corr(mU,T.predictions(TrainMaster,:));
    
    dC(:,:,tr) = differentialCore(Test_sC,Test_tC,mA,mB,cCore_S_composition_m,cCore_T_m,T.predictions(TestMaster,:));
end
%% 1) differntail random holdout with respect to genotytpe
for tr = 1:100
    T = f;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    UQ = unique(f.genoType);
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    
    
    
    %sC = [T.predictions];
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    Test_mU = Test_sC*mA;
    Test_mV = Test_tC*mB;
    UtoV = mU\mV;
    BETA = mA*UtoV*inv(mB);

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    cCore_S_composition_m = corr(mU,T.predictions(TrainMaster,:));
   
    
    dC(:,:,tr) = differentialCore(Test_sC,Test_tC,mA,mB,cCore_S_composition_m,cCore_T_m,T.predictions(TestMaster,:));
end
%% 2) differntail random holdout with respect to genotytpe
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];
UQ = unique(f.genoType);

for tr = 1:100
    T = M;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    tC = T.LAT;
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    sC = [T.predictions];
    
    
    %{
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    %}
    %[TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);    
    
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    
    if tr == 1
        BASEA = mA;
    else
        SGN(tr) = sign(BASEA(:,e)'*mA(:,e,tr));
        for e = 1:size(BASEA,2)
            if sign(BASEA(:,e)'*mA(:,e,tr)) == -1
                mA(:,e,tr) = -mA(:,e,tr);
                mB(:,e,tr) = -mB(:,e,tr);
                mU(:,e) = -mU(:,e);
                mV(:,e) = -mV(:,e);
            end
        end
        BASEA = mean(mA,3);
    end
    
    
    Test_mU = Test_sC*mA(:,:,tr);
    Test_mV = Test_tC*mB(:,:,tr);;
    UtoV = mU\mV;
    BETA = mA(:,:,tr)*UtoV*inv(mB(:,:,tr));

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    
    
    
    
    
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    
    
    
    
    
    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    %uK = mean(abs(CORE),3);
    %sK = std(abs(CORE),1,3)
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];
    
    for e = 1:2
        figure(tmp(e));
        bar(uK(e,:))
        hold on
        h = errorbar(uK(e,:),sK(e,:));
        set(h,'LineStyle','none');
        drawnow;
        hold off
    end
end
%% 3) geno typ ehold oout
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];
UQ = unique(f.genoType);
PRE = [];
ACC = [];
T = f;
for tr = 1:numel(UQ);
    
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    tC = T.LAT;
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    sC = [T.predictions];
    
    
    %{
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    %}
    %[TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);    
    TrainMaster = ~strcmp(T.genoType,UQ{tr});
    TestMaster = strcmp(T.genoType,UQ{tr});
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
    
    
    [gmA,gmB,gmr,gmU,gmV,gmstats] = canoncorr(sC,tC);
    gCORE = corr(gmU,T.predictions);
    
    if tr == 1
        BASEA = mA;
    else
        SGN(tr) = sign(BASEA(:,e)'*mA(:,e,tr));
        for e = 1:size(BASEA,2)
            if sign(BASEA(:,e)'*mA(:,e,tr)) == -1
                mA(:,e,tr) = -mA(:,e,tr);
                mB(:,e,tr) = -mB(:,e,tr);
                mU(:,e) = -mU(:,e);
                mV(:,e) = -mV(:,e);
            end
        end
        BASEA = mean(mA,3);
    end
    
    
    Test_mU = Train_sC*mA(:,:,tr);
    Test_mV = Test_tC*mB(:,:,tr);;
    UtoV = mU\mV;
    BETA = mA(:,:,tr)*UtoV*inv(mB(:,:,tr));

    
    
    PRE = [PRE;Test_sC*mA(:,:,tr)];
    ACC = [ACC;Test_tC*mB(:,:,tr)]
    
    
    
    
    
    
    
    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    
    
    
    
    
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    
    
    
    
    
    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    %uK = mean(abs(CORE),3);
    %sK = std(abs(CORE),1,3)
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];
    
    for e = 1:2
        figure(tmp(e));
        bar(uK(e,:))
        hold on
        h = errorbar(uK(e,:),sK(e,:));
        set(h,'LineStyle','none');
        drawnow;
        hold off
    end
end

%%
for e = 1:size(CORE,3)
    for cc = 1:size(CORE,1)
        if sign(squeeze(CORE(cc,:,e))*gCORE(cc,:)') == -1
            CORE(cc,:,e) = -CORE(cc,:,e);
        end
    end
end
K = abs(bsxfun(@minus,CORE,gCORE));
K = squeeze(K);
K = squeeze(sum(K,2));
sCORE = sum(gCORE,2);
K = bsxfun(@times,K,sCORE.^-1);
[s sidx] = sort(K(1,:));
sUQ = UQ(sidx);
figure
bar(s)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))
sortCORE = CORE(:,:,sidx);




for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    mv(u,:) = sum(gmU(fidx,:));
end
[s sidx] = sort(mv(:,1));
sUQ = UQ(sidx);
figure
bar(s)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))

figure;



%% JUNK








for e = 1:size(tC,2)
    CS = linspace(min(tC(TrainMaster,e)),max(tC(TrainMaster,e)),10);
    SIM = Core_T(e,:)'*CS;
    predict_tipAngle = PCA_BKPROJ(SIM',tE,tU);
    figure;plot(predict_tipAngle')
end


plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;


Prediction_tC = Test_sC*mA;
delta = Prediction_tC - Test_tC;
delta = sum(delta.*delta,2).^.5;

for e = 1:size(Prediction_tC,2)
    figure;
    plot(Prediction_tC(:,e),Test_tC(:,e),'.')
end

plot(Prediction_tC(:,1),Test_tC(:,1),'.')








for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end


% corr with predictions
[RHO PVAL] = corr(mU,T.predictions(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),T.tipAngle(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),T.tipAngle(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
 



%%
T = f;
close all

[mA,mB,mr,mU,mV,mstats] = canoncorr(T.predictions,T.shapeData);

[t_struct_C t_pval] = corr(mV,T.LAT);
for e = 1:size(mB,2)
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

[s_struct_C s_pval] = corr(mU,T.predictions);
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end
corr(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3));
figure
plot(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3),'.');
figure;
plot(mU(:,1),mV(:,1),'.');
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
ridx = find(any(abs(diff(f.tipAngle,1,2)) > 10*pi/180,2) | ...
            any(isnan(f.specData),2) | any(isnan(f.tipAngle),2) | ...
            any(abs(f.tipAngle(:,1)) > 45*pi/180,2) | ...
            any(abs(f.tipAngle(:,1) - f.tipAngle(:,end)) < 10*pi/180,2) | ...
            any(f.shapeData(:,4) > 3*10^5,2) | ...
            any(f.shapeData(:,5) > 400,2)  | ...
            any(abs(diff(f.length,1,2)) > 4,2)  | ...
            any(f.shapeDataL(:,end) > .5*10^9) | ...
            f.predictions(:,5) > f.predictions(:,6) | ... 
            isnan(f.xpos) ...
            );  
f.tipAngle(ridx,:) = [];
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
LAB.shapeOLD = fieldString(1005:1008);
LAB.shape = {'Height' 'Width' 'Depth' 'Volume'};
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
M.length = [f.length;Mf.length];
M.LAT = [f.LAT;Mf.LAT];
M.GR = [f.GR;Mf.GR];
M.shapeDataL = [f.shapeDataL;Mf.shapeDataL];
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



%% BEGIN ANALYSIS
%% A) no hold out for [composition] to [shape]
T = f;
close all
X = [T.predictions(:,1:end)];
Y = [T.shapeDataL(:,1:3)];

%{
X = zscore(X);
Y = zscore(Y);

[S X U Ex L ERR LAM] = PCA_FIT_FULL(X,size(X,2));
[S Y U Ey L ERR LAM] = PCA_FIT_FULL(Y,size(Y,2));
%}
%{
for e = 1:size(Ex,2)
    figure;
    bar(Ex(:,e));
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec));
    title(['Coff Spec:' num2str(e)]);
end

% PCA decomposition of seed shape
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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end


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
%% A) FIGURE 2 for paper
T = f;
close all
X = [T.specData];
Y = [T.shapeDataL(:,1:3)];

Y = zscore(Y);

[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);

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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end

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
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end
%% 

L = {'Weight' 'NMR' 'mg oil' 'oil percent'	'protein per' 'mg protein' 'total density' 'material density' 'total volume' 'material volume' 'air space'};
% biplot

c1 = [2^-.5 2^-.5]*[mU(:,1) mV(:,1)]';
c2 = [2^-.5 2^-.5]*[mU(:,2) mV(:,2)]';
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

%plot(mU(:,1)*norm(mA(:,1)),mU(:,2)*norm(mA(:,2)),'.');
plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');

[C_lessP] = PCA_REPROJ(lessP,Ex,Ux);
rm = [2];
L(rm) = [];
C_lessP(rm,:) = [];


% try to double check
for e = 1:size(mU,2)
    for l = 1:size(f.predictions,2)
        CHECK(e,l)= corr(f.predictions(:,l),mU(:,e));
    end
end
textMAG = 5;
C_lessP = CHECK;
hold all;
for e = 1:size(C_lessP,2)
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    %e1 = mA(:,1)/norm(mA(:,1));
    %p1 = e1'*C_lessP(e,:)';
    %e2 = mA(:,2)/norm(mA(:,2));
    %p2 = e2'*C_lessP(e,:)';
    %P = [p1 p2];    
    P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
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
        CHECK(e,l)= corr(Y(:,l),mU(:,e));
    end
end

L = {'Height' 'Width' 'Depth'};
C_lessP = CHECK;
for e = 1:size(C_lessP,2)    
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    %e1 = mA(:,1)/norm(mA(:,1));
    %p1 = e1'*C_lessP(e,:)';
    %e2 = mA(:,2)/norm(mA(:,2));
    %p2 = e2'*C_lessP(e,:)';
    %P = [p1 p2];    
    P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end

%{
CL = {'k.' 'c.' 'y.' 'm.' 'b.' 'r.' 'kd' 'cd' 'yd' 'md' 'bd' 'rd' ...
      'kp' 'cp' 'yp' 'mp' 'bp' 'rp' 'ks' 'cs' 'ys' 'ms' 'bs' 'rs' ...
      'kh' 'ch' 'yh' 'mh' 'bh' 'rh'};
UQ = unique(f.genoType);
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType));
    sub = mU(fidx,1:2);
    if size(sub,1) > 20
        [jS jC jU jE jL jERR jLAM] = PCA_FIT_FULL(sub,2);
        jLAM = diag(jLAM).^.5;
        ecc = axes2ecc(jLAM);
        [tx,ty] = ellipse1(mean(sub(:,1),1),mean(sub(:,2),1),[jLAM(1) ecc],180/pi*atan2(jE(2,1),jE(1,1)));
        plot(tx,ty,'k--','LineWidth',.2);
        plot(mean(sub(:,1)),mean(sub(:,2)),CL{u},'MarkerSize',7,'MarkerFaceColor',CL{u}(1));
    end
end
%}
%% A) FIGURE 2 for paper
T = f;
close all
X = [T.specData];
Y = [T.shapeDataL(:,1:3)];

Y = zscore(Y);

[S X Ux Ex L ERR LAM] = PCA_FIT_FULL(X,15);
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
    set(gca,'XTickLabel',LAB.shape,'XTick',1:numel(LAB.shape));
    title(['Coff Shape:' num2str(e)]);
end


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
    set(gca,'XTickLabel',LAB.spec,'XTick',1:numel(LAB.spec))
    title(['Core Spec:' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
    title(['Overall:' num2str(e)]);
end
%%
L = {'Weight' 'NMR' 'mg oil' 'oil percent'	'protein per' 'mg protein' 'total density' 'material density' 'total volume' 'material volume' 'air space'};
% biplot
c1 = [2^-.5 2^-.5]*[mU(:,1) mV(:,1)]';
c2 = [2^-.5 2^-.5]*[mU(:,2) mV(:,2)]';
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
% cv
plot(C(:,1),C(:,2),'.','MarkerSize',.2,'Color','k');



[C_lessP] = PCA_REPROJ(lessP,Ex,Ux);
rm = [2];
L(rm) = [];
C_lessP(rm,:) = [];



textMAG = 5;
MAG = 100;
hold all;
for e = 1:size(C_lessP,1)
    %C_lessP(e,:) = C_lessP(e,:)/norm(C_lessP(e,:));
    e1 = mA(:,1)/norm(mA(:,1));
    p1 = e1'*C_lessP(e,:)';
    e2 = mA(:,2)/norm(mA(:,2));
    p2 = e2'*C_lessP(e,:)';
    P = [p1 p2]/MAG;    
    %P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1),P(2),'LineWidth',3,'Color','k');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),LAB.spec{e},'Rotation',90+atan2(P(2),P(1))*180/pi,'FontSize',3);
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','k');
end
%%
MAG = 1;
C_lessP = eye(3);
L = {'Height' 'Width' 'Depth'};
for e = 1:size(C_lessP,1)  
    e1 = mB(:,1)/norm(mB(:,1));
    p1 = e1'*C_lessP(e,:)';
    e2 = mB(:,2)/norm(mB(:,2));
    p2 = e2'*C_lessP(e,:)';
    P = [p1 p2];    
    %P = C_lessP(1:2,e);
    %P = P/norm(P);        
    %P = P.*std(mU(:,1:2),1,1).^1;
    %quiver(0,0,p1,p2,'LineWidth',3);
    quiver(0,0,P(1)/MAG,P(2)/MAG,'LineWidth',3,'Color','r');
    P = P/norm(P);
    text(textMAG*P(1),textMAG*P(2),L{e},'Rotation',atan2(P(2),P(1))*180/pi,'FontSize',7,'Color','r');
    tX = linspace(0,textMAG*P(1),100);
    tY = linspace(0,textMAG*P(2),100);
    plot(tX,tY,'--','LineWidth',1,'Color','r');
end
















%% A) no hold out for [raw_spec] to [ICA_seedshape]
T = f;
close all
X = [T.predictions];
Y = [T.shapeDataL(:,1:3)];
Y = unmix';

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





%%
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];

PRE = [];
ACC = [];
T = M;
UQ = unique(T.genoType);
for tr = 1:numel(UQ);
    
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    tC = T.LAT;
    sC = [T.predictions];    
    
    TrainMaster = ~strcmp(T.genoType,UQ{tr});
    TestMaster = strcmp(T.genoType,UQ{tr});
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);    
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
    
    [gmA,gmB,gmr,gmU,gmV,gmstats] = canoncorr(sC,tC);
    gCORE = corr(gmU,T.predictions);
    
    
    PRE = [PRE;Test_sC*mA(:,:,tr)];
    ACC = [ACC;Test_tC*mB(:,:,tr)];
    
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];

end


for e = 1:size(CORE,3)
    for cc = 1:size(CORE,1)
        if sign(squeeze(CORE(cc,:,e))*gCORE(cc,:)') == -1
            CORE(cc,:,e) = -CORE(cc,:,e);
        end
    end
end
K = abs(bsxfun(@minus,CORE,gCORE));
K = squeeze(K);
K = squeeze(sum(K,2));
sCORE = sum(abs(gCORE),2);
K = bsxfun(@times,K,sCORE.^-1);
[sK sidx] = sort(K(1,:));
sUQ = UQ(sidx);
figure
bar(sK)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))
sortCORE = CORE(:,:,sidx);



% rank by 
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    mv(u,:) = sum(gmU(fidx,:));
end
[sMV sidx] = sort(mv(:,1));
sUQ = UQ(sidx);
figure
bar(sMV)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))

figure;hold all
for e = 1:size(mv,1)
    plot(abs(mv(e,1))/max(abs(mv(:,1))),abs(K(1,e))/max(abs(K(1,:))),'.');
end
legend(UQ)


%%
close all
plot(M.predictions(:,end-2),M.shapeDataL(:,end),'.')




%% try yet again
TEN = zeros(size(tC,2),size(sC,2),size(tC,1));
for e = 1:size(M.tipAngle,1)
    TEN(:,:,e) = tC(e,:)'*sC(e,:)*(sC(e,:)*sC(e,:)')^-1;
    e
end
%%
TENs = reshape(TEN,[size(TEN,1)*size(TEN,2) size(TEN,3)]);
[mS mC mU mE mL mERR mLAM] = PCA_FIT_FULL(TENs',3);
plot3(mC(:,1),mC(:,2),mC(:,3),'.')
%% 
sTENs = reshape(mS',size(TEN));
for e = 1:size(sC,1)
    tmpA = sTENs(:,:,e)*sC(e,:)';
    ang = PCA_BKPROJ(tmpA',tE,tU);
    plot(ang,'r');
    hold on;
    plot(M.tipAngle(e,:),'b');
    hold off
    drawnow
    pause(.5);
end

    










%% MASTER COMPLEX
T = M;
T = f;
T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st grouping - genotyp
T.groupingFactors(1).toGroup = 1;
T.groupingFactors(1).Groups = T.genoType;
T.groupingFactors(1).toDisplay = 1;
T.groupingFactors(1).Name = 'genoType';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd grouping - 50% whole
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
T.groupingFactors(2).toGroup = 0;
T.groupingFactors(2).Groups = holdOutGroups;
T.groupingFactors(2).toDisplay = 0;
T.groupingFactors(2).Name = 'random hold out across genotype';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3rd grouping - X percent % per genoType
UQ = unique(T.genoType);
perDrawTest = .75;
Train = zeros(size(T.specData,1),1);
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [sTrain, sTest] = crossvalind('HoldOut',size(fidx,1),perDrawTest);
    subI = find(sTrain);
    Train(fidx(subI)) = 1;
end
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
T.groupingFactors(3).toGroup = 0;
T.groupingFactors(3).Groups = holdOutGroups;
T.groupingFactors(3).toDisplay = 0;
T.groupingFactors(3).Name = 'random structured hold out accross genotype';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4th grouping kmeans
kidx = kmeans(T.specData,10);
for e = 1:numel(kidx)
    holdOutGroups{e} = num2str(kidx(e));
end
T.groupingFactors(4).toGroup = 0;
T.groupingFactors(4).Groups = holdOutGroups;
T.groupingFactors(4).toDisplay = 0;
T.groupingFactors(4).Name = 'kmeans display';


% configure pca for domain and codomain
type.PCAX.perform = 0;
type.PCAX.dim = 15;
type.PCAY.perform = 1;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
%type.learnMethod.which = 'cca';
type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
%type.learnMethod.which = 'lookup';
type.learnMethod.numComp = 15;
display.pca.perform.level = 0;
display.pca.perform.across = 1;
display.grouping = [];
display.grouping(1).grp = T.genoType;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO(T,holdOutGroups,type,display);


%% MASTER SIMPLE
T = f;
%T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));

perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
holdOutGroups = [T.genoType];


% configure pca for domain and codomain
type.PCAX.perform = 0;
type.PCAX.dim = 50;
type.PCAY.perform = 0;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
display.indivdual.perform = 0;
%type.learnMethod.which = 'cca';
%type.learnMethod.which = 'kcca';
type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
%type.learnMethod.which = 'lookup';
type.learnMethod.numComp = 15;
display.pca.perform = 1;
display.grouping = [];
display.grouping(1).grp = T.genoType;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO_OLD(T,holdOutGroups,type,display);

%% MASTER with loop over factors
cnt = 1;
RET = {};
for l = 2:1:30

    T = f;
    T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
    perDrawTest = .5;
    [Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
    for e = 1:numel(Train)
        holdOutGroups{e} = num2str(Train(e));
    end
    holdOutGroups = [T.genoType];

    % configure pca for domain and codomain
    type.PCAX.perform = 1;
    type.PCAX.dim = l;
    type.PCAY.perform = 1;
    type.PCAY.dim = 3;
    type.learnMethod.means = 0;
    type.learnMethod.which = 'cca';
    %type.learnMethod.which = 'pls';
    %type.learnMethod.which = 'net';
    %type.learnMethod.which = 'kMani';
    type.learnMethod.which = 'lookup';
    type.learnMethod.numComp = l;
    display.pca.perform = 0;
    display.indivdual.perform = 0;
    display.pca.domain = 'Spec Data';
    display.pca.codomain = 'Tip Angle';
    RET{cnt} = GO(T,holdOutGroups,type,display);
    cnt = cnt + 1;
    l
    
    
    RHO = [];
    for e = 1:numel(RET)    
        [rho pval] = corr(RET{e}.compPreGrpMeanY,RET{e}.compAccGrpMeanY);
        RHO(e,:) = diag(rho);
    end
    plot(RHO);
    drawnow
end
%% MASTER with shape
T = f;
T.tipAngle = bsxfun(@minus,T.tipAngle,T.tipAngle(:,1));
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(T.specData,1),perDrawTest);
for e = 1:numel(Train)
    holdOutGroups{e} = num2str(Train(e));
end
holdOutGroups = [T.genoType];
% configure pca for domain and codomain
type.PCAX.perform = 1;
type.PCAX.dim = 100;
type.PCAZ.perform = 0;
type.PCAZ.dim = 100;
type.PCAY.perform = 1;
type.PCAY.dim = 3;
type.learnMethod.means = 0;
type.learnMethod.which = 'cca';
%type.learnMethod.which = 'pls';
%type.learnMethod.which = 'net';
%type.learnMethod.which = 'kMani';
type.learnMethod.numComp = 60;
display.pca.perform = 1;
display.indivdual.perform = 1;
display.pca.domain = 'Spec Data';
display.pca.codomain = 'Tip Angle';
R = GO_SUM(T,holdOutGroups,type,display);
%% simple
X = M.specData;
X = bsxfun(@minus,X,mean(X,1));
X = bsxfun(@times,X,std(X,1,1));
Y = M.tipAngle;
Y = bsxfun(@minus,Y,mean(Y,1));
Y = bsxfun(@times,Y,std(Y,1,1));
sXY = (X'*Y);
for i = 1:size(X,2)
    for j = 1:size(Y,2)
        [RHO(i,j) PVAL(i,j)] = corr(X(:,i),Y(:,j));
    end
    i
end
%% explore again
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.tipAngle,3);
for e = 1:size(xS,1)
    plot(xS(e,:),'r')
    hold on
    plot(M.tipAngle(e,:),'b')
    hold off
    drawnow
    pause(.3)
end
%%
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.specData,3);
plot3(xC(:,1),xC(:,2),xC(:,3),'.','MarkerSize',3);
hold on
CL = {'r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','r','b','k','g','c','g','c','r','b','k','g','c'};
for u = 1:numel(UQ)
    idx = find(strcmp(M.genoType,UQ{u}));
    if numel(idx) > 3
        %plot3(xC(idx,1),xC(idx,2),xC(idx,3),[CL{u} 'o'],'MarkerSize',5);
        [xxS xxC xxU xxE xxL xxERR xxLAM] = PCA_FIT_FULL(xC(idx,:),3);
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,1),xxE(2,1),xxE(3,1),.1,'b');
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,2),xxE(2,2),xxE(3,2),.1,'g');
        quiver3(xxU(1),xxU(2),xxU(3),xxE(1,3),xxE(2,3),xxE(3,3),.1,'r');
    end
end
%% tray var model
for u = 1:numel(UQ
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(M.specData,3);






%% 1) spec --> shape -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-shape/';
mkdir(oPath);
Y = f.shapeData;
X = f.specData;
dimX = 10;
dimY = 9;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);


lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Top\_Area';
L{2} = 'Top\_Major';
L{3} = 'Top\_Minor';
L{4} = 'Front\_Area';
L{5} = 'Front\_Major';
L{6} = 'Front\_Minor';
L{7} = 'Side\_Area';
L{8} = 'Side\_Major';
L{9} = 'Side\_Minor';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end
%{
%% 2) spec --> tip --  create tip angle from spec data with hold out but group over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));
%X = gradient(X);
%X = gradient(X);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));
nG = 1;
dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
X = bsxfun(@minus,X,mean(X,1));
MV = figure;
timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(f.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

[wx,wy] = myCCA(xC(Train==1,:),yC(Train==1,:),3);
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),15);

% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([20],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}
for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};



        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];

oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%}



















%% 2) NAM ONLY - PLS - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY - CCA - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - PLS - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp(G,UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
%AUG = PCA_BKPROJ(BETA',xE,xU);
%inv(B)
%figure;plot(AUG');
%tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(X,1);
    dC = mean(X,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,X,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = uSPEC;
    %C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - CCA - spec --> tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end


for u = 1:numel(UQ)
    
    fidx = find(strcmp(G,UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end


%% 2) NAM ONLY PLS - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 3;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);

for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}));
    fidxTR = find(~strcmp([f.genoType],UQ{u}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};


        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY CCA - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = f.genoType;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);

for u = 1:numel(UQ)
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    [A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM CCA - spec --> tip --  create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = M.genoType;
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 50;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);



timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(G);

for u = 1:numel(UQ)
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    [A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end


%% 2) NAM ONLY - CCA - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = [f.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 5;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

% try clustering the X
kidx = kmeans(X,28);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    [A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    BETA = A*inv(B);
    
    fidx = find(strcmp([f.genoType],UQ{u}));
    fidxTR = find(~strcmp([f.genoType],UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) NAM ONLY - PLS - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
G = [f.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 400;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    %BETA = A*inv(B);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(uGCx,uGCy,15);
    BETA = BETA(2:end,:);
    
    fidx = find(strcmp(G,UQ{u}));
    %fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
    %% 5) spec --> MUT tip --  predict the classic muts from the NAM parents MEANS
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
xC = PCA_REPROJ(X,xE,xU);
yC = PCA_REPROJ(Y,yE,yU);
UQ = unique(Mf.genoType);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;



for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);


%% 2) CCM & NAM - CCA - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 5;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    [A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    BETA = A*inv(B);
    
    fidx = find(strcmp(G,UQ{u}));
    fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
figure;plot(AUG');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end
%% 2) CCM & NAM - PLS - spec --> tip --  MEANS -- create tip angle from spec data with hold out over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = M.tipAngle;
X = M.specData;
G = [M.genoType];
Y = bsxfun(@minus,Y,Y(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
dimX = 400;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
cl = [];
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;


for u = 1:numel(UQ)
    uGCx = [];
    uGCy = [];
    UQ = unique(G);    
    uidx = find(~strcmp(G,UQ{u}));
    genoTypeList = G(uidx);
    iUQ = unique(genoTypeList);
    for eu = 1:numel(iUQ)
        %fidx = find(strcmp(UQ{eu},genoTypeList));
        %uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        %uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
        fidx = find(strcmp(iUQ{eu},G));
        uGCx(eu,:) = mean(xC(fidx,:),1);
        uGCy(eu,:) = mean(yC(fidx,:),1);
        if any(strcmp(G(fidx),UQ{u}))
            STOP = 1
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(uGCx,uGCy);
    %BETA = A*inv(B);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(uGCx,uGCy,15);
    BETA = BETA(2:end,:);
    
    fidx = find(strcmp(G,UQ{u}));
    %fidxTR = find(~strcmp(G,UQ{u}));
    % cheat - for testing
    %fidxTR = [fidxTR;fidx];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clustered learning manifold
    %nG = 1;
    %lM = lManifold();
    %lM.setmodelCompX(dimX);
    %lM.setmodelCompY(dimY);
    %lM.addXY(xC(Train==1,:),yC(Train==1,:));
    %lM.setGroupN(nG);
    %lM.learn();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple canoncorr
    %[A,B,r,U,V,stats] = canoncorr(xC(fidxTR,:),yC(fidxTR,:));
    %BETA = A*inv(B);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple pls regress
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
    %BETA = BETA(2:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %net = feedforwardnet([15]);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % neural network
    %for e = 1:size(yC,2)
    %   net{e} = fitnet([10],'trainbr');
    %    net{e} = fitnet([20],'trainscg');
    %    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
    %end
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        %k = X(fidx,:)*BETA;
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        [RHO PVAL] = corr(k,yC(fidx,:));        
        fprintf(['Number of trials:' num2str(size(k,1)) '@' num2str(RHO(1)) ':' num2str(PVAL(1)) '\n']);                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        cl = [cl;u*ones(size(k,1),1)];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions - mean components
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];
        
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u1 = mean(k,1);
        s1 = std(k,1,1);
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % means for 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};

       
        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'...
      'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
    
    
    
    UQcl = unique(cl);
    for ui = 1:numel(UQcl)
        fidx = find(cl == ui);
        plot(subX(fidx),subY(fidx),CL{ui})
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions of means
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
AUG = PCA_BKPROJ(BETA',xE,xU);
inv(B)
figure;plot(AUG');
tot = [tot;AUG(3,:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    uSPEC = bsxfun(@plus,uSPEC,dC);
    C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
        axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end




%% 2) PLS - spec --> final_tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end
    %% 2) PLS - spec --> MUT final_tip --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);




Pre = [];
Acc = [];
uCPre = [];
uCAcc = [];


UQ = unique(Mf.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
%perDraw = .5;
%[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(Y(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


%% 2) PLS - spec --> avg_swingRate --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = f.tipAngle;
X = f.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);


timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
uCPre = [];
uCAcc = [];
scale = 180/pi;

UQ = unique(f.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);



fidxTR = find(Train==1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([f.genoType],UQ{u}) & Test==1);
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(yC(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end
    %% 2) PLS - spec --> MUT avg_swingRate --  create tip angle from spec data with hold out over rand
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Mf.tipAngle;
X = Mf.specData;
Y = Y(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prestage pca
%dimX = 10;
dimY = 1;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
xU = mean(X,1);
X = bsxfun(@minus,X,xU);




Pre = [];
Acc = [];
uCPre = [];
uCAcc = [];


UQ = unique(Mf.genoType);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build holdout index
%perDraw = .5;
%[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustered learning manifold
%nG = 1;
%lM = lManifold();
%lM.setmodelCompX(dimX);
%lM.setmodelCompY(dimY);
%lM.addXY(xC(Train==1,:),yC(Train==1,:));
%lM.setGroupN(nG);
%lM.learn();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple canoncorr
%[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%BETA = A*inv(B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple pls regress
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(fidxTR,:),yC(fidxTR,:),15);
%BETA = BETA(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% neural network
%for e = 1:size(yC,2)
%   net{e} = fitnet([10],'trainbr');
%    net{e} = fitnet([20],'trainscg');
%    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
%end




for u = 1:numel(UQ)
    
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % predict
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = X(fidx,:)*BETA;
        %k = X(fidx,:)*BETA(2:end,:);
        %k = net(xC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;Y(fidx,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather predictions
        uCPre = [uCPre;mean(k,1)];
        uCAcc = [uCAcc;mean(Y(fidx,:),1)];        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store genotype 
        genoVec{u} = UQ{u};
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at pca predictions 
for e = 1:size(Pre,2)
    subX = uCPre(:,e);
    subY = uCAcc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


%%
%%%%%%
%%


%%  2.25) spec --> final angle -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-finalTipAngle/';
mkdir(oPath);
Y = f.tipAngle(:,end);
X = f.specData;
dimX = 10;
dimY = 1;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Final Tip Angle';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end
%%  2.35) spec --> average swing rate -- nonlinear kmeans manifold with hold out with genotype as display
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-averageSwingRate/';
mkdir(oPath);
Y = mean(diff(f.tipAngle,1,2),2);
X = f.specData;
dimX = 10;
dimY = 1;
timer = 1;
scale = 180/pi;
nG = 1;
perDrawTest = .5;
[Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();
UQ = unique(f.genoType);
clear L
L{1} = 'Average Swing Rate';
A = [];
P = [];
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    if numel(uidx) >= 1
        subX = xC(uidx,:);
        subAY = Y(uidx,:);
        subY = lM.predict(subX);
        subY = PCA_BKPROJ(subY,yE,yU);
        P = [P;subY];
        A = [A;subAY];
    end
end

for e = 1:size(P,2)
    subX = P(:,e);
    subY = A(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    title([L{e} '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[oPath strrep(L{e},'\','') '.tif']);
end



%%      2.5) spec --> tip angle robust looped sample for predictions
Y = f.tipAngle;
Y = bsxfun(@minus,Y(:,1),Y);
X = f.specData;
dimX = 100;
dimY = 3;
scale = 180/pi;
nG = 1;
cor = figure;
pv = figure;
R = [];
pV = [];
for L = 1:50
    perDrawTest = .5;
    [Train, Test] = crossvalind('HoldOut',size(xC,1),perDrawTest);
    [xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
    [yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
    
    
    %net = feedforwardnet(20);
    %[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

    lM = lManifold();
    lM.setmodelCompX(dimX);
    lM.setmodelCompY(dimY);
    lM.addXY(xC(Train==1,:),yC(Train==1,:));
    lM.setGroupN(nG);
    lM.learn();
    
    %subY = net(xC(Test==1,:)')';
    subY = lM.predict(xC(Test==1,:));
    
    
    subY = PCA_BKPROJ(subY,yE,yU);
    subAY = Y(Test==1,:);
    [RHO,PVAL] = corr(subY,subAY);
    figure(cor);
    hold on
    plot(diag(RHO));
    R = [R diag(RHO)];
    figure(pv);
    hold on
    plot(diag(PVAL));
    pV = [pV diag(PVAL)];
end
figure(cor);
plot(mean(R,2),'k','LineWidth',5);
figure(pv);
plot(mean(pV,2),'k','LineWidth',5);

% display the groups in doamin
figure;
LEG = {};
for e = 1:size(lM.Ux,1)
    centerX = PCA_BKPROJ(lM.Ux(e,:),xE,xU);
    plot(centerX);
    LEG{e} = num2str(e);
    hold all
end
plot(mean(X),'k');
legend(LEG)
%p1 = anova1(lM.rawX(:,2),num2str(lM.kidx));
figure;
LEG = {};
for e = 1:size(lM.Uy,1)
    kidx = (lM.kidx == e);
    centerY = PCA_BKPROJ(lM.rawY(kidx,:),yE,yU);
        LEG{e} = num2str(e);
    ugY = mean(centerY,1);
    sgY = std(centerY,1,1)*size(centerY,1).^-.5;
    errorbar(ugY,sgY);
    hold all
end
legend(LEG);
%p1 = anova1(lM.rawY(:,1),num2str(lM.kidx));
%p2 = anova1(lM.rawY(:,2),num2str(lM.kidx));
%p3 = anova1(lM.rawY(:,3),num2str(lM.kidx));
plot(mean(Y),'k');
%% 3) shape --> tip --  create tip angle from spec data with hold out but group over genotype

close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.shapeData;
Y = bsxfun(@minus,Y(:,1),Y);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 9;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);


timer = .1;
YP = [];
XP = [];
Pre = []
Acc = [];
UQ = unique(f.genoType);
scale = -180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% neural network
net = feedforwardnet(3);
net.trainParam.max_fail = 20;
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        k = lM.predict(xC(fidx,:));
        k = net(xC(fidx,:)')';
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];




        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
    end

end





for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['SHAPE --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 3) shape X spec --> tip --  create tip angle from spec data with hold out but group over genotype

close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-tip/';
mkdir(outPath);
disp = 1;


Y = f.tipAngle;
X = f.shapeData;
Z = f.specData;
Y = bsxfun(@minus,Y(:,1),Y);


% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

dimX = 8;
dimY = 3;
dimZ = 8;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[zS zC zU zE zL zERR zLAM] = PCA_FIT_FULL(Z,dimZ);
X = [xC zC];
dimX = dimX+dimZ;
%[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
xC = X;


net = network();
net.numInputs = size(xC,2);
net.numLayers = 3;


% neural network
net = feedforwardnet(9);
net = fitnet(9);
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');


timer = .1;
YP = [];
XP = [];
Pre = []
Acc = [];
UQ = unique(f.genoType);
scale = -180/pi;


% clustered learning manifold
nG = 1;
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();




for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        k = lM.predict(xC(fidx,:));
        k = net(xC(fidx,:)')'; 
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];




        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
    end

end





for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['SHAPE --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end


[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 4) shape --> tip -- nonlinear kmeans manifold with hold out on genotype
Y = f.tipAngle;
X = f.shapeData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 10;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
PX = [];
PY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end
%% 5) spec --> MUT tip --  predict the classic muts from the NAM parents
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;
%p = randperm(size(X,1));
%Y = Y(p,:);
Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 10;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
myC = PCA_REPROJ(Mf.tipAngle,yE,yU);
mxC = PCA_REPROJ(Mf.specData,xE,xU);

timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% neural network
%net = feedforwardnet([15]);
net = fitnet([20]);
[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        k = lM.predict(mxC(fidx,:));
        
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;myC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 5) PLS spec --> MUT tip --  predict the classic muts from the NAM parents
%{
tmp = f;
f = Mf;
Mf = tmp;
%}
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;


Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 100;
dimY = 3;
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
myC = PCA_REPROJ(Mf.tipAngle,yE,yU);

% simple pls regress
X = bsxfun(@minus,X,mean(X,1));
mxC = bsxfun(@minus,Mf.specData,mean(f.specData,1));
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),5);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);


timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;


%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %k = lM.predict(mxC(fidx,:));
        k = mxC(fidx,:)*BETA;
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        %if strcmp(UQ{u}(end),'t')
            % gather predictions
            Pre = [Pre;k];
            Acc = [Acc;myC(fidx,:)];
        %end
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);
        
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
              %{  
                figure(MV);
                errorbar(u2*scale,s2*scale,'b');
                hold on
                errorbar(u1*scale,s1*scale,'m');
            %}

                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                hold off
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                %if mod(u,2) == 0
                %    hold off
                %end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% 5) PLS spec --> MUT tip --  predict the classic muts from the NAM parents
close all
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/shape-MUTtip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
MY = Mf.tipAngle;
p = randperm(size(X,1));
%Y = Y(p,:);
Y = bsxfun(@minus,Y,Y(:,1));
MY = bsxfun(@minus,MY,MY(:,1));
% build holdout index
perDraw = .01;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);

nG = 1;
dimX = 100;
dimY = 3;
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);

% simple pls regress
X = bsxfun(@minus,X,mean(X,1));
mxC = bsxfun(@minus,Mf.specData,mean(f.specData,1));
UX = [];
UY = [];
UQ = unique(f.genoType);
for u = 1:numel(UQ) 
    fidx = find(strcmp([f.genoType],UQ{u}));
    UX = [UX;mean(X(fidx,:),1)];
    UY = [UY;mean(yC(fidx,:),1)];
end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),10);
%[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(UX,UY,10);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);


timer = .02;
YP = [];
XP = [];
Pre = []
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(Mf.genoType);
scale = 180/pi;


%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([10],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}


for u = 1:numel(UQ)
    fidx = find(strcmp([Mf.genoType],UQ{u}));
    
    if ~isempty(fidx)
        %k = lM.predict(mxC(fidx,:));
        k = mxC(fidx,:)*BETA;
        %k = net(mxC(fidx,:)')';
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        if strcmp(UQ{u}(end),'t')
            % gather predictions
            Pre = [Pre;k];
            Acc = [Acc;myC(fidx,:)];
        end
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(MY(fidx,:),1);
        s2 = std(MY(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};

        if disp
            try
                if mod(u,2) == 0
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'b');
                    hold on
                    errorbar(u1*scale,s1*scale,'m');
                end
                
                if mod(u,2) == 1
                    figure(MV);
                    errorbar(u2*scale,s2*scale,'k');
                    hold on
                    errorbar(u1*scale,s1*scale,'r');
                end
                
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                if mod(u,2) == 0
                    hold off
                end
                waitforbuttonpress();
                if ~isempty(outPath)
                    fn = [outPath strrep(UQ{u},filesep,'') '.tif'];
                    saveas(gca,fn);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end

%figure;
%plotregression(Acc,Pre)

[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);





%% try nonlinear kmeans manifold with hold out on genotype - predict means and ind from means
Y = f.tipAngle;
X = f.specData;
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;

for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(subX);

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(3);
        
        lMStore{u} = lM;
    end
end

[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end


%% 2) PLS spec --> tip --  create tip angle from spec data with hold out but group over genotype
outPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/spec-tip/';
mkdir(outPath);
disp = 1;

Y = f.tipAngle;
X = f.specData;
Y = bsxfun(@minus,Y,Y(:,1));
%X = gradient(X);
%X = gradient(X);

% build holdout index
perDraw = .5;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
%Train = ones(size(Train));
nG = 1;
dimX = 100;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
X = bsxfun(@minus,X,mean(X,1));
MV = figure;
timer = .1;
YP = [];
XP = [];
Pre = [];
Acc = [];
uPre = [];
uAcc = [];
UQ = unique(f.genoType);
scale = 180/pi;

% clustered learning manifold
lM = lManifold();
lM.setmodelCompX(dimX);
lM.setmodelCompY(dimY);
lM.addXY(xC(Train==1,:),yC(Train==1,:));
lM.setGroupN(nG);
lM.learn();

% simple canoncorr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
BETA = A*inv(B);


% simple pls regress
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X(Train==1,:),yC(Train==1,:),15);
xC = X;
xU = 0;
xE = 1;
BETA = BETA(2:end,:);

% neural network
%net = feedforwardnet([15]);
%[net tr] = train(net,xC(Train==1,:)',yC(Train==1,:)');

%{
for e = 1:size(yC,2)
    %net{e} = fitnet([10],'trainbr');
    net{e} = fitnet([20],'trainscg');
    [net{e} tr] = train(net{e},xC(Train==1,:)',yC(Train==1,e)');
end
%}
for u = 1:numel(UQ)
    fidx = find(strcmp([f.genoType],UQ{u}) & Test == 1);
    
    if ~isempty(fidx)
        %k = lM.predict(xC(fidx,:));
        % predict with test
        %k = xC(fidx,:)*A*inv(B);
        k = xC(fidx,:)*BETA;
        %k = X(fidx,:)*BETA;
        %k = net(xC(fidx,:)')';
        
        %{
        k = [];
        for i = 1:3
            k = [k net{i}(xC(fidx,:)')'];
        end
        %}
        
        
        % gather predictions
        Pre = [Pre;k];
        Acc = [Acc;yC(fidx,:)];
        
        % back project
        k = PCA_BKPROJ(k,yE,yU);

        u1 = mean(k,1);
        s1 = std(k,1,1);
        
        u2 = mean(Y(fidx,:),1);
        s2 = std(Y(fidx,:),1,1);

        YP = [YP;k];
        XP = [XP;Y(fidx,:)];
    
        % 
        uPre = [uPre;u1];
        uAcc = [uAcc;u2];
        
        
        % store genotype 
        genoVec{u} = UQ{u};



        if disp
            try
                figure(MV);
                errorbar(u2*scale,s2*scale,'k');
                hold on
                errorbar(u1*scale,s1*scale,'r');
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(timer);
                hold off
                %waitforbuttonpress();
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
                
            catch

            end
        end
    end

end

for e = 1:size(Pre,2)
    subX = Pre(:,e);
    subY = Acc(:,e);
    [RHO PVAL] = corr(subX,subY);
    MX = max(([subX(:);subY(:)]));
    MN = min(([subX(:);subY(:)]));
    figure;
    plot(subX,subY,'.');
    axis([MN MX MN MX]);
    hold on
    plot(linspace(MN,MX,3),linspace(MN,MX,3),'r');
    ylabel('Predicted');
    xlabel('Actual');
    title(['Spec --> Tip Angle -nGroups' num2str(nG) ' - PC' num2str(e) '--R value--' num2str(RHO) '--pVal--' num2str(PVAL)]);
    saveas(gca,[outPath 'PC-' num2str(e) '.tif']);
end



% pearson corr
[RHO,PVAL] = corr(XP,YP);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);

% augment the prediction vectors
%AUG = PCA_BKPROJ(lM.corVecA{1}',lM.Ex(:,:,1),lM.Ux(1,:));
%AUG = PCA_BKPROJ(BETA',xE,xU);
AUG = BETA';
figure;plot(AUG');
oPath = '/mnt/spaldingdata/nate/communications/labmeetings/jeff_nate/LATtipAngleChange/';
mkdir(oPath);
% look at oil on tip angle
for e = 1:size(AUG,1)
    uSPEC = mean(f.specData,1);
    dC = mean(f.specData,1);
    uSPEC = uSPEC - dC;
    mag = 1;
    TEST = AUG(e,:)*bsxfun(@minus,f.specData,dC)';
    oil = AUG(e,:);
    delta = mag*std(TEST)*oil*(oil*oil')^-1;
    delta = delta'*[-1 0 1];
    %delta = delta'*[0];
    uSPEC = bsxfun(@plus,uSPEC,delta');
    %uSPEC = bsxfun(@plus,uSPEC,dC);
    %C = PCA_REPROJ(uSPEC,xE,xU);
    %C = squeeze(lM.pVec(1,:,:))'*C';
    C = uSPEC;
    C = BETA'*C';
    C = PCA_BKPROJ(C',yE,yU);
    figure;

    plot(C'*scale);
    axis([0 61 -30 90]);
    nVar = ['Lat' num2str(e)];
    title(['Change in Tip Angle wrt ' nVar]);
    xlabel('Time (fr)');
    ylabel('Angle (degees)');
    saveas(gca,[oPath nVar '.tif']);
end



%% pie graphs
for g = 1:nG
    sidx = find(Train==1);
    fidx = find(lM.kidx == g);
    gs = f.genoType(sidx(fidx));
    for u = 1:numel(UQ)
        tG(u) = sum(strcmp(gs,UQ{u}));
    end
    tG = tG / sum(tG);
    figure;
    pie(tG);
        title(['Group' num2str(g)]);
    legend(UQ);
end
%% try nonlinear kmeans manifold with hold out on genotype
Y = f.tipAngle;
X = f.specData;
%Y = bsxfun(@minus,Y,Y(:,1));
%Y = diff(Y,1,2);
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
UQ = unique(f.genoType);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        lM.addXY(xC(uidx,:),yC(uidx,:));
        lM.setGroupN(nG);
        lM.learn();
    
        %net = fitnet(20);
        %net = train(net,xC(uidx,:)',yC(uidx,:)');
        
        subX = xC(Test,:);
        subAY = Y(Test,:);
        
        subY = lM.predict(subX);
        %subY = net(subX')';
        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        [RHO,PVAL] = corr(subY,subAY);
        figure;plotyy(1:size(RHO,1),diag(RHO),1:size(RHO,1),diag(PVAL));
        title(UQ{u});
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(.1);
    end
end
[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX',uPY');
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    for u = 1:numel(UQ)
        plot(mean(cPX(pUQ==u,e)),mean(cPY(pUQ==u,e)),'r*');
    end
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title([num2str(RHO) '--' num2str(PVAL)]);
end



[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);
%% BELOW predicitons are with genotype means for building model with hold out
%% spec --> tip -- nonlinear kmeans manifold with hold out on genotype spec
Y = f.tipAngle;
X = f.specData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
title('R value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(RHO)) max(diag(RHO))])
figure;plot(diag(PVAL));
title('p Value');
hold on;
plot(diag(PVAL) < .05,'r');
axis([0 61 min(diag(PVAL)) max(diag(PVAL))]);


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end

% p
for g = 1:size(lM.pVec,1)
    pV = PCA_BKPROJ(squeeze(lM.pVec(g,:,:))',xE,xU);
    figure;
    plot(pV')
    title([num2str(g)]);
end
% look at 
for g = 1:numel(lM.corVecA)
    pV = PCA_BKPROJ(squeeze(lM.corVecA{g})',xE,xU);
    figure;
    plot(pV')
    title([num2str(g)]);
end
%{
hold on
plot(pC(1,:),'k');
pC(1,:)*lM.corVecA{g}'
%}
% 
for g = 1:numel(lM.corVecB)
    pV = PCA_BKPROJ(squeeze(lM.corVecB{g})',yE,yU);
    figure;
    plot(pV')
    title([num2str(g) '-- with means added' ]);
end
%
for g = 1:numel(lM.corVecB)
    pV = PCA_BKPROJ(squeeze(lM.corVecB{g})',yE,yU);
    pV = bsxfun(@minus,pV,yU);
    figure;
    plot(pV')
    title([num2str(g) '-- with means subtracted' ]);
end
%% spec x shape -- > tip -- nonlinear kmeans manifold with hold out on genotype 
Y = f.tipAngle;
X = f.shapeData;
Z = f.specData;
%Y = bsxfun(@minus,Y(:,1),Y);
dimX = 3;
dimY = 3;
dimZ = 8;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
[zS zC zU zE zL zERR zLAM] = PCA_FIT_FULL(Z,dimZ);
X = [xC zC];
dimX = dimX+dimZ;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;
for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(mean(subX,1));

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;mean(yC(Test,:),1)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(1);
    end
end
% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end
%% try nonlinear kmeans manifold with hold out on genotype - predict means and ind from means
Y = f.tipAngle;
X = f.specData;
dimX = 8;
dimY = 3;
nG = 1;
scale = 180/pi;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);
PX = [];
PY = [];
uPX = [];
uPY = [];
cPX = [];
cPY = [];
pUQ = [];
uGCx = [];
uGCy = [];
MI = figure;

for u = 1:numel(UQ)-1
    uidx = find(~strcmp([f.genoType],UQ{u}));
    
    genoTypeList = f.genoType(uidx);
    UQ = unique(genoTypeList);
    for eu = 1:numel(UQ)
        fidx = find(strcmp(UQ{eu},genoTypeList));
        uGCx(eu,:) = mean(xC(uidx(fidx),:),1);
        uGCy(eu,:) = mean(yC(uidx(fidx),:),1);
    end
    
    Test = find(strcmp([f.genoType],UQ{u}));
    if numel(Test) >= 10
        %uidx = [Test;uidx];
        lM = lManifold();
        lM.setmodelCompX(dimX);
        lM.setmodelCompY(dimY);
        
        lM.addXY(uGCx,uGCy);
        lM.setGroupN(nG);
        lM.learn();
    
        subX = xC(Test,:);
        subAY = Y(Test,:);
        subY = lM.predict(subX);

        pUQ = [pUQ;u*ones(numel(Test),1)];
        cPX = [cPX;yC(Test,:)];
        cPY = [cPY;subY];
        
        subY = PCA_BKPROJ(subY,yE,yU);

        uP = mean(subY,1);
        sP = std(subY,1,1);
        uA = mean(subAY,1);
        sA = std(subAY,1,1);
        
        PX = [PX;subAY];
        PY = [PY;subY];
        
        uPX = [uPX;uA];
        uPY = [uPY;uP];
        
        
        figure;
        plot(subY'*scale,'b--');
        hold on
        plot(subAY'*scale,'g--');
        errorbar(uP*scale,sP*scale,'r','LineWidth',2);
        hold on
        errorbar(uA*scale,sA*scale,'k','LineWidth',2);
        hold off
        axis([0 61 -30 90]);
        %axis([0 61 -2 4]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        
        figure(MI);
        errorbar(uP*scale,sP*scale,'r');
        hold on
        errorbar(uA*scale,sA*scale,'k');
        hold off
        axis([0 61 -30 90]);
        title(UQ{u});
        xlabel('Time (fr)');
        ylabel('Tip Angle (deg)');
        drawnow
        pause(3);
        
        lMStore{u} = lM;
    end
end

[RHO,PVAL] = corr(PX,PY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));

[RHO,PVAL] = corr(uPX,uPY);
figure;plot(diag(RHO));
figure;plot(diag(PVAL));


% create PC predictions plots
for e = 1:size(cPX,2)
    figure;
    plot(cPX(:,e),cPY(:,e),'.');
    hold on
    plot(linspace(-100,100),linspace(-100,100),'r');
    AX = 5*std(cPX(:,e));
    axis([-AX AX -AX AX]);
    [RHO,PVAL] = corr(cPX(:,e),cPY(:,e));
    title(['PC' num2str(e) '--R value--' num2str(RHO) '--p value--' num2str(PVAL)]);
end


%% look at indivduals tip angles per genotype
for u = 1:numel(UQ)
    uidx = find(strcmp([f.genoType],UQ{u}));
    figure;
    plot(f.tipAngle(uidx,:)','g--');
    hold on
    plot(mean(f.tipAngle(uidx,:)),'k','LineWidth',5)
end
%% look at oil,protien,density,volume
X = f.tipAngle;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
for p = 1:size(f.predictions,2)
    for e = 1:size(xC,2)
        figure
        plot(f.predictions(:,p),xC(:,e),'.');
        title([num2str(p) '--' num2str(e)]);
    end
end
%% create prediction plot(s) of values for tip angle vs oil,protien, density, volume
% build holdout index
perDraw = .01;
Y = f.tipAngle;
X = f.predictions;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
plotSet = Train;
dimX = 4;
dimY = 3;
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,dimX);
[tS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,dimY);    


%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
% predict with test
k = xC(plotSet==1,:)*A*inv(B);
for i = 1:size(yC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% kmeans on oil and plot tip angle
NG = 5;
oil_idx = kmeans(f.predictions(:,2),NG);
clear LEG
for g = 1:NG
    uidx = find(oil_idx==g);
    uT = mean(f.tipAngle(uidx,:),1);
    sT = std(f.tipAngle(uidx,:),1,1)*numel(uidx)^-.5;
    errorbar(uT,sT);
    hold all
    uP(g) = mean(f.predictions(uidx,2));
    LEG{g} = num2str(uP(g));
end
legend(LEG)




%% normalize
for e = 1:size(pC,1)
    npC(e,:) = pC(e,:)/norm(pC(e,:));
end
npC*npC'



















%% 0) load shape data
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
S = schemas(conn);
T = tables(conn);
P = columns(conn);
sid = 'kernel_3d';
for e = 1:size(P,1)
    if strcmp(P{e},sid)
        cNames = P{e,2};
        sidNUM = e;
    end
end
%% 0.1) construct query for kernel 3D
q = ['select * from kernel_3d '...    
    'join kernels on kernel_3d.kernel_id = kernels.id ' ...
    'join kernel_plates on kernels.plate_id = kernel_plates.id '...
    'join population_lines on kernel_plates.population_line_id = population_lines.id'];
results = fetch(conn,q);
%% 0.2) put kernel features into hashmap
import java.util.HashMap;
kernelVec = HashMap();
ROWnum = 6;
COLnum = 8;
ROWvec = {'A','B','C','D','E','F'};
COLvec = {'1','2','3','4','5','6','7','8'};
for e = 1:size(results,1)
    
    try
        NUMpos = results{e,45}-1;
        
        
        if NUMpos > ROWnum*COLnum
            NUMpos = NUMpos - ROWnum*COLnum;
        end
        plateN = results{e,50};
       
        if strcmp(plateN(end),'A') | strcmp(plateN(end),'B')  
            plateN = [plateN(1:end-1) '-' plateN(end)];    
        end
        
        
        
        colN = mod(NUMpos,COLnum) + 1;
        colV = ['_' COLvec{colN}];

        rowN = floor(NUMpos/COLnum) + 1;
        rowV = ['_' ROWvec{rowN}];
        
        
    catch ME

        
    end
    
    
    wellN = [rowV colV];

    key1 = [plateN '*' wellN];
    if strfind(plateN,'07S-MO001')
        plateN
    end
    %fprintf([key1 '-->' wellN '-->' num2str(NUMpos) '\n']);
    kernelVec.put(key1,results(e,3:42));
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 1) load csv masterlist
[masterList, result]= readtext('/mnt/spaldingdata/florida/files/nir_master_list.csv');
%% 2) load oil strach protein prdict vectors
[pC] = readtext('/mnt/spaldingdata/nate/COEFF.csv');
pC = cell2mat(pC);
pC(:,end) = [];
for i = 1:size(pC,1)
    npC(i,:) = pC(i,:)/norm(pC(i,:));
end
%% 3) make platename -->genotype key value store from masterlist
import java.util.HashMap;
K1 = 7;
V1 = 11;
P2G = HashMap();
for e = 1:size(masterList,1)
    key = masterList{e,K1};
    value = masterList{e,V1};    
    if ~isempty(value)
        P2G.put(key,value);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 4) dig for tip angle data
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction/loganSPOOL/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% 6) open tip angles
cnt = 1;
DT = [];
WNT = {};
PLT = {};
for i = 1:numel(FileList)
    try
        [p n ext] = fileparts(FileList{i});    
        fidx = strfind(n,'_');

        plateName = n(fidx(1)+1:fidx(2)-1);
        % factor
        wellNames = n(fidx(2)+1:end);
        D = csvread(FileList{i});
        wellNames = ['_' wellNames '_'];
        fidx = strfind(wellNames,'_');
    
        if numel(fidx)-1 == size(D',2)
            DT = [DT D'];
            for i = 1:numel(fidx)-1
                WNT{cnt} = wellNames(fidx(i)+1:fidx(i+1)-1);
                PLT{cnt} = plateName;
                cnt = cnt + 1;
            end
        end
    catch ME
        ME
      
    end
end
%% 7) clean tip angles
% clean via derivative
ridx = find(any(abs(diff(DT,1,1)) > 20/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
% clean via total bend
ridx = find(any(abs(DT(1,:) - DT(end,:)) < 20/180*pi,1));
DT(:,ridx) = [];
PLT(ridx) = [];
WNT(ridx) = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 8) run import init from mongoDB
import phytoG.locked.Bpersist.Bos.implementations.*
import com.mongodb.*;
import java.util.Map.*;
import java.util.HashMap;
import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
import phytoG.locked.BdataObjects.BbioObjects.maize.spectraData;
oStore = OStore_mdb();
oStore.accessResource();
oStore.setCollection('specData');
%% 9) import spec and pull kernel shape data from hashmap
f.specData = [];
f.tipAngle = [];
f.kernelData = [];
f.weightData = [];
f.genoType = {};
f.wellNumber = {};
f.plateName = {};
for e = 1:numel(WNT)
    qMap = HashMap();
    pn = PLT{e};
    qMap.put('_k_m._pnode._k_m._plateName',pn);
    wn = ['_' WNT{e}(1) '_' WNT{e}(2)];
    qMap.put('_k_m._pnode._k_m._wellName',wn);
    cursor = oStore.search(qMap);
    itr = cursor.iterator();
    if itr.size()
        n = itr.next();
        n = spectraData(n);
        
        spec = n.getSpectrum();
        
        
        
        %%%%%%%
        % pull spec data from array
        specData = [];
        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        
        %%%%%%%
        % get kernel shape data
        sKey = [pn '*' wn];
        kVec = kernelVec.get(sKey);
        kernelData = [];
        if ~isempty(kVec)
            for i = 1:kVec.size()
                kernelData(i,1) = kVec(i);
            end
        end
        
        
        next_specData = specData;
        next_kernelData = kernelData;
        next_weight = n.getWeight();
        next_tipAngle =  DT(:,e);
        next_genoType = P2G.get(pn(1:end-2));
        next_wellName = wn;
        next_plateName = pn;
        
        %%%%% check for valid entries
        if ~isempty(next_specData) && ...
           ~isempty(next_kernelData) && ...
           ~isempty(next_tipAngle) && ...
           ~isempty(next_genoType) && ...
           ~isempty(next_wellName) && ...
           ~isempty(next_plateName) && ...
           ~isempty(next_weight)
           
           
            f.specData = [f.specData next_specData];
            f.tipAngle = [f.tipAngle next_tipAngle];
            f.weightData = [f.weightData str2num(char(next_weight))];
            f.kernelData = [f.kernelData next_kernelData];
            f.genoType{end+1} = next_genoType;
            f.wellNumber{end+1} = next_wellName;
            f.plateName{end+1} = next_plateName;
            fprintf(['plate: ' pn ' well: ' wn ' found! \n']);
        else
            fprintf(['plate: ' pn ' well: ' wn ' incomplete Data found! \n']);   
        end
        
        
        
    else
        fprintf(['plate: ' pn ' well: ' wn ' NOT found! \n']);
    end
end
f.specData(1:3,:) = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 10) get subset of kernel feature data - area and lengths
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
%% 10.5) CHOICE - look at histogram(s)
for i = 1:size(f.kernelSub,1)
    figure
    hist(f.kernelSub(i,:))
end
%% 11) clean data
% clean on area
ridx = find(f.kernelSub(1,:) > 4*10^5);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);

ridx = find(f.kernelSub(2,:) > 500);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
%% 12) get subsubset of kernel data - the lengths
sIDX = [11 30 31] - 2;
f.kernelLengths = f.kernelData(sIDX,:);
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% JUNK
disp = 1;
phytoK = {};
ridx = [];
for e = 1:numel(f.plateName)
    try
        topName = translateToKernelImage(f.plateName{e},f.wellNumber{e});
        kernelTopExtract2(topName,disp);
    end
end
%% XX-13) look at image via lookups
disp = 1;
phytoK = {};
ridx = [];
for e = 1:numel(f.plateName)
    try
        topName = translateToKernelImage(f.plateName{e},f.wellNumber{e});
        [height width area contour] = kernelTopExtract2(topName,disp);
        phytoShapes = [height width area];
        k = [f.plateName{e} '*' f.wellNumber{e}];
        phytoK(e).key = k;
        phytoK(e).value = phytoShapes;
        phytoK(e).curve = contour;
        fprintf([k '\n']);
        ridx(e) = 0;
    catch ME
        phytoK{e} = [];
        ridx(e) = 1;
        fprintf(['error@' topName '\n']);
    end
end
%% JUNK
ridx = find(ridx);
f.specData(:,ridx) = [];
f.tipAngle(:,ridx) = [];
f.weightData(:,ridx) = [];
f.kernelData(:,ridx) = [];
f.genoType(ridx) = [];
f.wellNumber(ridx) = [];
f.plateName(ridx) = [];
% recreate
sIDX = [9 10 11 19 20 21 29 30 31] -2;
f.kernelSub = f.kernelData(sIDX,:);
sIDX = [11 30 31] - 2;
f.kernelLengths = f.kernelData(sIDX,:);
%% XX-13) do inserts for phytoKernel
f.phytoKernel = [];
for e = 1:numel(f.plateName)
     k = [f.plateName{e} '*' f.wellNumber{e}];
     fidx = find(strcmp({phytoK.key},k));
     f.phytoKernel(e,:) = phytoK(e).value;
end
f.phytoKernel = f.phytoKernel';
%% BREAK
%% look at lengths by lengths
close all
for i = 1:size(f.kernelLengths,1)
    for j = 1:size(f.kernelLengths,1)
        figure;
        plot(f.kernelLengths(i,:),f.kernelLengths(j,:),'.')
        title([num2str(i) '--' num2str(j)]);
    end
end
%% 11.2) create weights from lengths with hold out
close all
disp = 1;
X = f.kernelLengths';
Y = f.weightData';
perDraw = .7;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);

[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,3);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,1);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
[mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),1);
% predict with test
k = sC(Test==1,:)*A*inv(B);
mk = sC(Test==1,:)*mA*inv(mB);


k = PCA_BKPROJ(k,tE,tU);
mk = PCA_BKPROJ(mk,tE,tU);


plot(k,Y(Test==1),'.')
figure;
plot(mk,Y(Test==1),'.');

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%% MUST) flip tip angle
f.tipAngle = -f.tipAngle;





%% 10) hold out
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',6);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end
%%%%%%%%
%% BREAK
%% BREAK
%% BREAK
%% lets look at volume predictions
close all
offset = 0;
plot3(f.kernelSub(offset+1,:),f.kernelSub(offset+2,:),f.kernelSub(offset+3,:),'.');
figure;
plot(f.kernelSub(offset+1,:),f.kernelSub(offset+2,:).*f.kernelSub(offset+3,:),'.');
figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.kernelSub(offset+1:offset+3,:)',3);
plot3(sC(:,1),sC(:,2),sC(:,3),'.');
%% 10.5) look at corr between weights and area,lengths and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelSub';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};

for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
end
legend(UQ);
%% 10.5) look at corr between weights and length and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelLengths';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'along length of cob - front minor','kernel height - top major','along circ of cob - top minor'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Lengths');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 125;
%% 10.5) look at corr between weights and phytoLengths and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.phytoKernel';

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'along length of cob - front minor','kernel height - top major','along circ of cob - top minor'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Lengths');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 125;
%% 10.5) look at corr between weights and area and legend on genotype
close all
disp = 1;
Y = f.weightData';
X = f.kernelData(27,:)';  % only look at top area

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
LABELS = {'top area'};
for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
    xlabel('Kernel Area');
    ylabel('Kernel Weights');
    title([LABELS{i} ' - vs weight'])
end
% mini script for lookup
fidx1 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) < 200;
fidx2 = strcmp([f.genoType],'MS71') & f.kernelLengths(2,:) > 200;
fidx3 = strcmp([f.genoType],'MS71') & f.kernelLengths(1,:) < 125;
%% 10.6) look at volumes and density vs weight(s)
% volume(S)
Y = f.weightData';
X = [f.kernelData(8,:).*f.kernelData(9,:).*.5.*(f.kernelData(18,:)+f.kernelData(28,:));...
     f.kernelData(18,:).*f.kernelData(19,:).*.5.*(f.kernelData(8,:)+f.kernelData(29,:));...
     f.kernelData(28,:).*f.kernelData(29,:).*.5.*(f.kernelData(19,:)+f.kernelData(9,:))]';
plot3(X(:,1),X(:,2),X(:,3),'r.');
title('Test each volumes vs themselves');
axis equal
%X = mean(X,2).^(1/3);
% calc density from volume and mass
X = mean(X,2);
X = [X Y.*X.^-1];

CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};

for i = 1:size(X,2)
    figure;
    hold on
    UQ = unique([f.genoType]);
    for u = 1:numel(UQ)
        fidx = strcmp([f.genoType],UQ{u});
        plot(X(fidx,i),Y(fidx),CL{u});
    end
    legend(UQ);
end
%% BREAK
%% 11) create tip angle from spec data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
for L = 1:1
    % loop over number of basis vectors
    for num = 15

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %{
        %%%%%%%%%%%%%%%
        % pls regression        
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from spec data AND shape data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = [f.specData;f.kernelSub]';
for L = 1:1
    for num = 15

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}

        %{
        %%%%%%%%%%%%%%%
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from kernel-Lengths with NO hold out
%close all
disp = 1;
Y = f.tipAngle';
X = [f.kernelLengths]';
for L = 1:1
    for num = 3

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        % predict
        k = sC*A*inv(B);

        
        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}

        %{
        %%%%%%%%%%%%%%%
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);
            
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            
        end
        
    end
end
%% 11) create tip angle from spec data AND shape data with NO hold out
close all
disp = 1;
Y = f.tipAngle';
X = [f.specData]';
%X2 = [f.specData;f.kernelSub]';
X2 = [f.specData]';
for L = 1:1
    for num = 5

        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [xS2 sC2 sU2 sE2 sL2 sERR2 sLAM2] = PCA_FIT_FULL(X2,num);
        [xS2s sC2s sU2s sE2s sL2s sERR2s sLAM2s] = PCA_FIT_FULL(f.kernelSub',6);        
        sC2 = [sC2 sC2s];
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

        %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(:,yi),tmpX);
            k = [k tmpX*b];
        end
        %}


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        [A2,B2,r2,U2,V2,stats2] = canoncorr(sC2,tC);
        % predict
        k = sC*A*inv(B);
        k2 = sC2*A2*inv(B2);


        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}


        k = PCA_BKPROJ(k,tE,tU);
        k2 = PCA_BKPROJ(k2,tE,tU);
        
        UQ = unique([f.genoType]);

        for u = 1:numel(UQ)
            % find the genotype
            fidx = find(strcmp([f.genoType],UQ{u}));
            u1 = mean(k(fidx,:));
            u12 = mean(k2(fidx,:));
            u2 = mean(Y(fidx,:),1);


            s1 = std(k(fidx,:),1,1);
            s12 = std(k2(fidx,:),1,1);
            s2 = std(Y(fidx,:),1,1);

            dif(u) = mean(u1*u2');
            if disp
                try
                    errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    hold on
                    errorbar(u12*-180/pi,s12*-180/pi,'g');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                    axis([0 61 -30 90])
                    title(UQ{u});
                    pause(1);
                    hold off
                catch

                end
            end
            %}
        end
        DIF(num) = mean(dif);
    end
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.2) create tip angle from spec data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_spec/';
disp = 1;
Y = f.tipAngle';
X = f.specData';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,15);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.2) create tip angle from shape data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_shape/';
disp = 1;
Y = f.tipAngle';
X = f.kernelLengths';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,3);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),3);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.2) create tip angle from shape data with hold out
close all
outPath = '/mnt/spaldingdata/nate/communications/meetings/NSF_talk_Edgar_2013/predict_tip_angle_from_shapeANDspec/';
disp = 1;
Y = f.tipAngle';
X = f.kernelLengths';
Z = f.specData';
perDraw = .5;
% build holdout index
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
for L = 1:1

    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL([X Z],7);    
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

    
    
    %%%%%%%%%%%%%%%
    % perform corr
    [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
    [mA,mB,mpx,mpy] = myCCA(sC(Train==1,:),tC(Train==1,:),3);
    % predict with test
    k = sC(Test==1,:)*A*inv(B);
    mk = sC(Test==1,:)*mA*inv(mB);


    k = PCA_BKPROJ(k,tE,tU);
    mk = PCA_BKPROJ(mk,tE,tU);
    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)

        fidx = find(strcmp([f.genoType(Test==1)],UQ{u}));
        wfidx = find(strcmp([f.genoType],UQ{u}));


        u1 = mean(k(fidx,:));
        mu1 = mean(mk(fidx,:));
        u2 = mean(Y(wfidx,:),1);


        s1 = std(k(fidx,:),1,1);
        ms1 = std(mk(fidx,:),1,1);
        s2 = std(Y(wfidx,:),1,1);


        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');
                hold on
                %errorbar(mu1*-180/pi,ms1*-180/pi,'b');
                errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                xlabel('Time (fr)');
                ylabel('Tip Angle (deg)');
                pause(1);
                hold off
                if ~isempty(outPath)
                    saveas(gca,[outPath UQ{u} '.tif']);
                end
            catch

            end
        end
   
    end
end
%% 11.3) create hold out via genotype
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
UQ = unique([f.genoType]);
num = 15;
for L = 1:3
    for u = 1:numel(UQ)

        Test = strcmp([f.genoType],UQ{u});
        Train = ~ Test;


        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);    

        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train,:),tC(Train,:));
        % predict with test
        k = sC(find(Test==1),:)*A*inv(B);


        k = PCA_BKPROJ(k,tE,tU);

        u1 = mean(k);
        u2 = mean(Y(Test==1,:),1);

        s1 = std(k,1,1);
        s2 = std(Y(Test==1,:),1,1);

        
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                axis([0 61 -30 90])
                title(UQ{u});
                pause(1);
                hold off
            catch

            end
        end
    end
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.7) create prediction plot(s) of values for tip angle vs specdata
% build holdout index
perDraw = .01;
Y = f.tipAngle';
X = f.kernelSub';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

num = 6;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
%[pA,pB] = myPLS1(sC(Train==1,:),tC(Train==1,:));
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for tip angle vs shapedata
% build holdout index
perDraw = .01;
Y = f.tipAngle';
X = f.specData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for tip angle vs kernel Lengths
% build holdout index
close all
perDraw = .01;
Y = f.tipAngle';
X = f.kernelLengths';
numX = 3;
numY = 3;
plotSet = Train;
[Train, Test] = crossvalind('HoldOut', size(X,1),perDraw);
[xS xC xU xE xL xERR xLAM] = PCA_FIT_FULL(X,numX);
[yS yC yU yE yL yERR yLAM] = PCA_FIT_FULL(Y,numY);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(xC(Train==1,:),yC(Train==1,:));
%[A,B,U,V] = myCCA(xC(Train==1,:),yC(Train==1,:),5);
% predict with test
k = xC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = xC(plotSet==1,i);
    %yP = k(:,i);
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo');
    %axis equal
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    %plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    %axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for spec vs shapedata
% build holdout index
perDraw = .01;
Y = f.kernelSub';
X = f.specData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

numX = 6;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,numX);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% 11.7) create prediction plot(s) of values for weight vs shapedata
% build holdout index
perDraw = .01;
numX = 6;
numY = 1;
X = f.kernelSub';
Y = f.weightData';
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
plotSet = Train;

[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,numX);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,numY);

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
%[A,B,U,V] = myCCA(sC(Train==1,:),tC(Train==1,:),5);
% predict with test
k = sC(plotSet==1,:)*A*inv(B);

for i = 1:size(tC,2)
    xP = U(:,i);
    yP = V(:,i);
    
    %xP = sC(plotSet==1,i);
    %yP = k(:,i);
    
    
    [RHO,PVAL] = corr(xP,yP);
    
    figure;
    hold on;
    plot(xP,yP,'mo')
    axis equal
    
    
    mpre = min(xP);
    mC = min(yP);    
    wMIN = min(mpre,mC);
    
    Mpre = max(xP);
    MC = max(yP);    
    wMAX = max(Mpre,MC);
    
    
    %title([num2str(r(i)) '-' num2str(stats.p(i))]);
    title([num2str(r(i))  '----' num2str(RHO) '-----' num2str(PVAL) '-----' num2str(stats.p(i))]);
    hold on
    plot(linspace(wMIN,wMAX,100),linspace(wMIN,wMAX,100),'r');
    axis([wMIN wMAX wMIN wMAX]);
end
%% BREAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 11.8) create oil, starch, protien from tip angles via CCA
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;    
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via CCA by genotype
% build holdout index
f.ospw = pC*f.specData;
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,num);    

%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);    
    figure;
    hold on
    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        realV = mean(f.ospw(i,fidx));
        preV = mean(k(fidx,i));        
        plot(realV,preV,'.');        
    end
    title(TITLE{i})
    axis equal
end
%% 11.8) create oil, starch, protien from tip angles via pls
% build holdout index
perDraw = .3;
f.ospw = pC*f.specData;
[Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
num = 11;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    

% pls regression
[XL,YL,XS,YS,BETA] = plsregress(tC,sC,5);
% predict
k = [ones(size(tC,1),1) tC]*BETA;


k = PCA_BKPROJ(k,sE,sU);
k = k*pC';
TITLE = {'oil','starch','protien','weight'};
for i = 1:size(f.ospw,1);
    figure;        
    %plot(f.ospw(i,Test==1),k(:,i),'.')
    plot(f.ospw(i,:),k(:,i),'.')
    title(TITLE{i})
    axis equal
end
%% 11.9) create oil, starch, protien from tip angles directly
% build holdout index
f.ospw = pC*f.specData;
X = f.ospw';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,4);


%%%%%%%%%%%%%%%
% perform corr
[A,B,r,U,V,stats] = canoncorr(sC,tC);
% predict with test
k = tC*B*inv(A);
k = PCA_BKPROJ(k,sE,sU);

TITLE = {'oil','starch','protien','weight'};
UQ = unique([f.genoType]);
for i = 1:size(f.ospw,1);        
    figure;
    hold on
    
    realV = (f.ospw(i,:));
    preV = (k(:,i));        
    plot(realV,preV,'.');        

    title(TITLE{i})
    axis equal
end
%% test myCCA
Y = f.tipAngle';
X = f.specData';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);
[A,B,r,U,V,stats] = canoncorr(sC,tC);
[wx,wy] = myCCA(sC,tC,3);
%%

close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
fitX = f.specData';
fitY = f.tipAngle';
cnt = 1;
for num = 3:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        
        
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(fitX,num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(fitY,3);
        
        
    %{
        %%%%%%%%%%%%%%%
        % linear regression
        k = [];
        tmpX = [ones(size(sC,1),1) sC];
        for yi = 1:size(tC,2)
            [b,bint,r,rint,stats] = regress(tC(Train==1,yi),tmpX(Train==1,:));
            k = [k tmpX*b];
        end
        %}
        


        %%%%%%%%%%%%%%%
        % perform corr
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));
        % predict
        %predictY = sC*A*inv(B);
        predictY = sC*A*B';



        %{
        % pls regression
        [XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
        % predict
        k = [ones(size(sC,1),1) sC]*BETA;
        %}
        
        predictY = PCA_BKPROJ(predictY,tE,tU);
        
        
        trTOP = sum(predictY(Train==1,:).*fitY(Train==1,:),2);
        trBOTX = sum(predictY(Train==1,:).*predictY(Train==1,:),2);
        trBOTY = sum(fitY(Train==1,:).*fitY(Train==1,:),2);
        trCORR = trTOP.*trBOTX.^-.5.*trBOTY.^-.5;
        
        tsTOP = sum(predictY(Test==1,:).*fitY(Test==1,:),2);
        tsBOTX = sum(predictY(Test==1,:).*predictY(Test==1,:),2);
        tsBOTY = sum(fitY(Test==1,:).*fitY(Test==1,:),2);
        tsCORR = tsTOP.*tsBOTX.^-.5.*tsBOTY.^-.5;
        
        
        trainR(loop) = mean(trCORR);
        testR(loop) = mean(tsCORR);        
        
    end
    
    trU(cnt) = mean(trainR);
    trS(cnt) = std(trainR);
    
    tsU(cnt) = mean(testR);
    txS(cnt) = std(testR);
    
    cnt = cnt + 1;
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end
%% try correlation via genotype
close all
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',5);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);


UQ = unique([f.genoType]);
tipFig = figure;
specFig = figure;
ccFig = figure;
ccFigT = figure;
ccFigS = figure;
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'}
for u = 1:numel(UQ)
    -180/pi
    
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData(:,fidx)',8);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle(:,fidx)',8);
    
    
    % get the tip angle
    tmpTip = f.tipAngle(:,fidx)';
    % get the spec data
    tmpSpec = f.specData(:,fidx)';
    % get the means
    tU = mean(tmpTip,1);
    sU = mean(tmpSpec,1);
    % get the std
    tS = std(tmpTip,1,1)*size(tmpTip,1)^-.5;
    xS = std(tmpSpec,1,1)*size(tmpSpec,1)^-.5;
    
    if numel(fidx) > 1
        %[A,B,r,U,V,stats] = canoncorr(sC(fidx,:),tC(fidx,:));
        [A,B,r,U,V,stats] = canoncorr(sC,tC);
        
        
        
        
        figure(ccFigS);
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        plot(sData);
        hold all;
        
        figure(ccFigT);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(tData);
        hold all;
        
        
        LEG{u} = [UQ{u} '--' num2str(r(1)) '--' num2str(numel(fidx)) '--' num2str(stats.p(1))];
        figure(ccFig);
        scatter(U(:,1),V(:,1),CL{u});axis equal
        hold all
    else
        figure(ccFigT);
        scatter(0,0,CL{u});
        LEG{u} = [UQ{u} '---NA'];
    end
    
    
    figure(tipFig);
    errorbar(tU,tS);
    hold all
    
    figure(specFig);
    errorbar(sU,xS);
    hold all
    
    
end
figure(tipFig);
legend(UQ);
figure(specFig);
legend(UQ);
figure(ccFig);
legend(LEG);
figure(ccFigS);
legend(UQ);
figure(ccFigT);
legend(UQ);
%% try genotype cluster

f.ospw = pC*f.specData;


figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(sC(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,stats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);

scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
f.ospw = pC*f.specData;
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
%X = f.ospw';
for num = 15
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(Y(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(Y(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}%% oil starch protien

f.ospw = pC*f.specData;

figure;
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',4);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

UQ = unique([f.genoType]);
CL = {'r.','b.','g.','k.','m.','c.','y.',...
      'r*','b*','g*','k*','m*','c*','y*',...
      'ro','bo','go','ko','mo','co','yo',...
      'r^','b^','g^','k^','m^','c^','y^'};
utC = [];
usC = [];
for u = 1:numel(UQ)
    % find the uth group
    fidx = find(strcmp([f.genoType],UQ{u}));
    if numel(fidx) > 2
        usC = [usC;mean(f.ospw(fidx,:))];
        utC = [utC;mean(tC(fidx,:))];
    end
end


[A,B,r,U,V,st
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
endats] = canoncorr(usC,utC);
sData = PCA_BKPROJ(A(:,1)',sE,sU);
tData = PCA_BKPROJ(B(:,1)',tE,tU);
figure;plot(sData);

scatter(U(:,1),V(:,1));axis equal;
axis([-3 3 -3 3]);
%% please try
f.ospw = pC*f.specData;
close all
disp = 1;
Y = f.tipAngle';
X = f.specData';
%X = f.ospw';
for num = 15
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(X,num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,5);    
    
    
    [A,B,r,U,V,stats] = canoncorr(sC,tC);
    k = sC*A*inv(B);
    k = PCA_BKPROJ(k,tE,tU);



    UQ = unique([f.genoType]);

    for u = 1:numel(UQ)
        fidx = find(strcmp([f.genoType],UQ{u}));
        u1 = mean(k(fidx,:));
        u2 = mean(Y(fidx,:),1);

        
        s1 = std(k(fidx,:),1,1);
        s2 = std(Y(fidx,:),1,1);
        
        dif(u) = mean(u1*u2');
        if disp
            try
                errorbar(u1*-180/pi,s1*-180/pi,'r');hold on;errorbar(u2*-180/pi,s2*-180/pi,'k')
                %axis([0 61 -30 90])
                pause(1);
                hold off
            catch

            end
        end
        %}
    end
    DIF(num) = mean(dif);
end
sK = PCA_BKPROJ(A',sE,sU);
%{
for e = 1:size(k,1)
    plot(k(e,:),'k');
    hold on
    plot(f.tipAngle(:,e),'r');
    hold off
    drawnow
    pause(.4)
end
%}
%% size matters
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',size(f.specData',2));
newData = (diag(diag(sLAM).^-.5)*sC')';
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(newData,5);
%%
ng = 1;
kidx = kmeans(f.tipAngle',ng);
kidx = kmeans(f.specData',ng);
%% init first tip angle to zero
f.tipAngle2 = bsxfun(@minus,f.tipAngle,f.tipAngle(1,:));
%f.tipAngle2 = diff(f.tipAngle,1,1);
%% 
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',8);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',8);
T = [sC tC];
ng = 1;
kidx = kmeans(T,ng);C = PCA_REPROJ(,Ey,Uy);
%% try upbulk 2 term
for i = 1:size(sC,2)
    for j = i:size(sC,2)
        sC = [sC sC(:,i).*sC(:,j)];
    end
end
sE%% explore nonlinear
plot3(sC(:,1),sC(:,2),tC(:,2),'.')
%% pls regression
[xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',11);
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
[XL,YL,XS,YS,BETA] = plsregress(sC,tC,3);
Y0 = [ones(size(sC,1),1) sC]*BETA;
for e = 1:size(Y0,1)
    plot(Y0(e,:),'r');
    sData = PCA_BKPROJ(A(:,1)',sE,sU);
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% clustered cannon corr --> 
CL = {'r.','b.','g.','k.','m.'};
figure;
for num = 10
    [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);
    for e = 1:ng
        [A,B,r,U,V,stats] = canoncorr(sC(kidx==e,:),tC(kidx==e,:));
        sData = PCA_BKPROJ(A(:,1)',sE,sU);
        tData = PCA_BKPROJ(B(:,1)',tE,tU);
        plot(U(:,1),V(:,1),CL{e});axis equal
        hold on;
        R(num) = r(1);
    end
    hold off;
end

axis([-4 4 -4 4])
figure;plot(sData);
figure;plot(tData);
figure;plot(-f.tipAngle);
figure;plot(f.specData);
%%
close all
perDraw = 1/3;
testVec = 1;
trU = [];
tsU = [];
trS = [];
txS = [];
for num = 1:70
    
    trainR = [];
    testR = [];
    
    for loop = 1:30
        [xS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData',num);
        [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle',3);

        % build holdout index
        [Train, Test] = crossvalind('HoldOut', size(f.specData,2),perDraw);
        % build mapping
        [A,B,r,U,V,stats] = canoncorr(sC(Train==1,:),tC(Train==1,:));

        
        x = sC(Train==1,:)*A(:,testVec);
        y = tC(Train==1,:)*B(:,testVec);
        trainR(loop) = mean(x.*y);

        x = sC(Test==1,:)*A(:,testVec);
        y = tC(Test==1,:)*B(:,testVec);
        testR(loop) = mean(x.*y);
    end

    
    trU(num) = mean(trainR);
    trS(num) = std(trainR);
    
    tsU(num) = mean(testR);
    txS(num) = std(testR);
    
    errorbar(trU,trS,'k');
    hold on
    errorbar(tsU,txS,'r');    
    drawnow
    
end


%%
pth = '/mnt/spaldingdata/Bessie/forPoster/';
csvwrite([pth 'tipAngle.csv'],f.tipAngle);
%%
csvwrite([pth 'corr.csv'],[U(:,1),V(:,1)]);
%% view via eye
OP = tData;
mea = OP*f.tipAngle;
[J sidx] = sort(mea);

hold on;
for e = 1:numel(sidx)
    plot(tS','b');
    hold on
    plot(tS(sidx(e),:),'r');
    hold off
    drawnow
end
%% predict spec
A = sC'/tC';
sim = A*tC';
M = PCA_BKPROJ(sim',sE,sU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.specData(:,e),'b');
    hold off
    pause(.3)
end
%% predict tip angle
A = tC'/sC';
sim = A*sC';
M = PCA_BKPROJ(sim',tE,tU);
for e = 1:size(M,1)
    plot(M(e,:),'r');
    hold on
    plot(f.tipAngle(:,e),'b');
    hold off
    pause(.3)
end
    %%
    cnt = 1;
    while (itr.hasNext())
        n = itr.next();
        n = spectraData(n);
        spec = n.getSpectrum();

        for i = 1:spec.size()
            specData(i,1) = str2num(spec.get(i-1));
        end
        DS(:,cnt) = specData;
        WNS{cnt} = n.getPlateName();
        PNS{cnt} = n.getWellName();
        cnt = cnt + 1
    end
%% joy again
[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(f.specData,tC,30);
%[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(xsC,xtC,3);
Y = [ones(size(tC,1),1) f.specData]*BETA(1:end,:);
%Y = [ones(size(tC,1),1) xsC]*BETA(1:end,:);
plot(tC(:,1),Y(:,1),'.');
%plot(XS(:,1),YS(:,1),'.')
%%
  [A,B,r,U,V,stats] = canoncorr(f.specData,f.tipAngle);

  [A,B,r,U,V,stats] = canoncorr(specSTACK,tipSTACK);
    %[A,B,r,U,V,stats] = canoncorr(xC,yC);

figure;
plot(U(:,1),V(:,1),'.');
title([num2str(r(1)) '--' num2str(stats.p(1))]);
        % predict with test
        k = subX(Test==1,:)*A*inv(B);
%%
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData,20);
UQ = unique(f.genoType);
tipSTACK = [];
specSTACK = [];
G = {};
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType)==1);
    if numel(fidx) > 8
        [utS utC utU utE utL utERR utLAM] = PCA_FIT_FULL(tC(fidx,:),3);
        for bv = 1:size(utE,2)
            utE(:,bv) = utE(:,bv)*sign(utE(bv,bv));
            diag(utE);
        end
        [utC] = PCA_REPROJ(tC(fidx,:),utE,utU);
        utC = utC*diag(utL.^-1);
        %tipSTACK = [tipSTACK; utC*diag(utL.^-1)];
        %tipSTACK = [tipSTACK; utC];
        tipSTACK = [tipSTACK; tC(fidx,:)];
        
        
        [usS usC usU usE usL usERR usLAM] = PCA_FIT_FULL(sC(fidx,:),3);
        for bv = 1:size(usE,2)
            usE(:,bv) = usE(:,bv)*sign(usE(bv,bv));
        end
        [usC] = PCA_REPROJ(sC(fidx,:),usE,usU);
        usC = usC*diag(usL.^-1);
        %specSTACK = [specSTACK; usC*diag(usL.^-1)];
        %specSTACK = [specSTACK; usC];
        specSTACK = [specSTACK; sC(fidx,:)];
        G = [G ;f.genoType(fidx)];
    end
end
%%
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(f.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(f.specData,15);
UQ = unique(f.genoType);
tipSTACK = [];
specSTACK = [];
G = {};
for u = 1:numel(UQ)
    fidx = find(strcmp(UQ{u},f.genoType)==1);
    if numel(fidx) > 6
        %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(sC(fidx,:),tC(fidx,:),3);
        [A,B,r,U,V,stats] = canoncorr(sC(Train,:),tC(Train,:));
        xsC(fidx,:) = U;
        xtC(fidx,:) = V;
        G = UQ{u};
    end
end
specSTACK = sxC;
tipSTACK = stC;

%%
UQ = unique(G);
DSe(1) = figure;
DSe(2) = figure;
DSe(3) = figure;

DSr(1) = figure;
DSr(2) = figure;
DSr(3) = figure;
ux = [];
uy = [];
uRx = [];
uRy = [];
[mA,mB,mr,mU,mV,mstats] = canoncorr(specSTACK,tipSTACK);

for u = 1:numel(UQ)
    Test = find(strcmp(UQ{u},G));
    Train = find(~strcmp(UQ{u},G));
    
    sUQ = setdiff(G,UQ{u});

    gX = [];
    gY = [];
    for u2 = 1:numel(sUQ)
        idx2 = find(strcmp(sUQ{u2},G));
        gX(u2,:) = mean(specSTACK(idx2,:),1);
        gY(u2,:) = mean(tipSTACK(idx2,:),1);
    end
    %Test = find(strcmp(UQ{u},f.genoType));
    %Train = find(~strcmp(UQ{u},f.genoType));
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(specSTACK(Train,:),tipSTACK(Train,:),3);
    %A = XL;
    %B = YL;
    %[XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(f.specData(Train,:),tC(Train,:),3);
    
    %plot(XS(:,1),YS(:,1),'.')
    
    [A,B,r,U,V,stats] = canoncorr(specSTACK(Train,:),tipSTACK(Train,:));
    MA(:,:,u) = A;
    MB(:,:,u) = B;
    %[A,B,r,U,V,stats] = canoncorr(gX,gY);
    u1 = mean(specSTACK(Train,:),1);
    u2 = mean(tipSTACK(Train,:),1);
    %BETA = A*inv(B);
    
    %ygit = [ones(numel(Test),1) specSTACK(Test,:)]*BETA;
    %ygit = [ones(numel(Test),1) f.specData(Test,:)]*BETA;
    %ygit = [specSTACK(Test,:)]*BETA;
    xfit = bsxfun(@minus,specSTACK(Test,:),u1)*A;
    yfit = bsxfun(@minus,tipSTACK(Test,:),u2)*B;    
    ux = [ux ;mean(xfit,1)];
    uy = [uy ;mean(yfit,1)];
    
    for e = 1:size(yfit,2)
        figure(DSe(e))
        hold on
        %plot(ygit(:,e),tC(Test,e),'b.');
        plot(xfit(:,e),yfit(:,e),'b.');        
    end
    %ygit = [specSTACK(Train,:)]*BETA(2:end,:);
    %ygit = [ones(numel(Train),1) specSTACK(Train,:)]*BETA;
    %ygit = [ones(numel(Train),1) f.specData(Train,:)]*BETA;
    %ygit = [specSTACK(Train,:)]*BETA;
    xfit = bsxfun(@minus,specSTACK(Train,:),u1)*A;
    yfit = bsxfun(@minus,tipSTACK(Train,:),u2)*B;    
    uRx = [uRx ;mean(xfit,1)];
    uRy = [uRy ;mean(yfit,1)];
    for e = 1:size(yfit,2)
        figure(DSr(e))
        hold on
        %plot(ygit(:,e),tC(Test,e),'b.');
        plot(xfit(:,e),yfit(:,e),'r.');
        plot(ux(end,e),uy(end,e),'*');
    end
end
for e = 1:size(yfit,2)
    figure(DSr(e))
    hold on
    plot(ux(:,e),uy(:,e),'k*');
    plot(uRx(:,e),uRy(:,e),'g*');
end
[RHO,PVAL] = corr(ux,uy);
%% master of joy - hold out 50% genotype of each one
T = f;
initT = T.tipAngle(:,1);
T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
UQ = unique(T.genoType);
perDraw = .5;
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,15);
%[shS shC shU shE shL shERR shLAM] = PCA_FIT_FULL(T.shapeData,3);
%sC = [sC shC];
TrainMaster = [];
TestMaster = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
    TrainMaster = [TrainMaster;Train];
    TestMaster = [TestMaster;Test];
end
TrainMaster = find(TrainMaster);
TestMaster = find(TestMaster);
u_sC = mean(sC(TrainMaster,:),1);
Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
u_tC = mean(tC(TrainMaster,:),1);
Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
[mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
Test_mU = Test_sC*mA;
Test_mV = Test_tC*mB;
BETA = mA*inv(mB);
plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end

% corr with predictions
[RHO PVAL] = corr(mU,T.predictions(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),T.tipAngle(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),T.tipAngle(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
%% joy master with mutants - predict from above
T = Mf;
initT = T.tipAngle(:,1);
T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
tC_mut = PCA_REPROJ(T.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(T.specData,sE,sU);
UQ = unique(T.genoType);
pre = figure; 
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
   
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    figure;
    errorbar(pUTA,pSTA,'r');
    hold on
    axis([0 61 0 2])
    errorbar(UTA,STA,'k');
    axis([0 61 0 2])
    hold off
    title(UQ{u});
end
%% scatter plots for mutants
initT = Mf.tipAngle(:,1);
Mf.tipAngle = bsxfun(@minus,Mf.tipAngle,initT);
tC_mut = PCA_REPROJ(Mf.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(Mf.specData,sE,sU);
UQ = unique(Mf.genoType);
preC(1) = figure; 
preC(2) = figure; 
preC(3) = figure; 
close all
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
tmpX = [];
tmpY = [];
UtmpX = [];
UtmpY = [];
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(mean(tC_mut(fidxw,e)),mean(predict(:,e)),CL{u},'MarkerSize',7);
       
       plot(linspace(-2,2,100),linspace(-2,2,100),'r');
    end
    UtmpX = [UtmpX;mean(tC_mut(fidxw,:),1)];
    UtmpY = [UtmpY;mean(predict,1)];
    
    tmpX = [tmpX;tC_mut(fidxw,:)];
    tmpY = [tmpY;predict];
    
end
%{
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predict = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*BETA;
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(tC_mut(fidxw,e),predict(:,e),[CL{u}(1) '.']);
       plot(linspace(-10,10,100),linspace(-10,10,100),'r');
    end
end
%}
%% scatter plots for mutants
initT = Mf.tipAngle(:,1);
Mf.tipAngle = bsxfun(@minus,Mf.tipAngle,initT);
tC_mut = PCA_REPROJ(Mf.tipAngle,tE,tU);
sC_mut = PCA_REPROJ(Mf.specData,sE,sU);
UQ = unique(Mf.genoType);
preC(1) = figure; 
preC(2) = figure; 
preC(3) = figure; 
close all
CL = {'r*' 'm*' 'b*' 'g*' 'k*' 'y*' 'c*' 'r^' 'm^' 'b^' 'g^' 'k^' 'y^' 'c^' ...
      'ro' 'mo' 'bo' 'go' 'ko' 'yo' 'co' 'r.' 'm.' 'b.' 'g.' 'k.' 'y.' 'c.'};
tmpX = [];
tmpY = [];
UtmpX = [];
UtmpY = [];
for u = 1:numel(UQ)
    fidxw = find(strcmp(Mf.genoType,UQ{u}));    
    predictSPEC = bsxfun(@minus,sC_mut(fidxw,:),u_sC)*mA;
    predictTIP = bsxfun(@minus,tC_mut(fidxw,:),u_tC)*mB;
    
    
    for e = 1:size(predict,2)
       figure(preC(e));
       hold on
       plot(mean(predictSPEC(:,e)),mean(predictTIP(:,e)),CL{u},'MarkerSize',7);
       
       plot(linspace(-2,2,100),linspace(-2,2,100),'r');
    end
    
    UtmpX = [UtmpX;mean(predictSPEC,1)];
    UtmpY = [UtmpY;mean(predictTIP,1)];
    
    tmpX = [tmpX;predictSPEC];
    tmpY = [tmpY;predictTIP];
    
end
%% look at model
cv_tipAngle = PCA_BKPROJ(mB',tE,tU);
plot(cv_tipAngle(1:2,:)');
cv_specData = PCA_BKPROJ(mA',sE,sU);
non_lat_spec = predictVectors*bsxfun(@minus,cv_specData,dX)';
for e = 1:size(non_lat_spec,2)
    figure;
    bar(non_lat_spec(:,e));
end
%% try tip angle or length whole from only the predictions
T = f;%M;
%initT = T.tipAngle(:,1);
%T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
%UQ = unique(T.genoType);
perDraw = .05;
%Y = gradient(T.length);
%Y = T.length;
[tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(Y,3);
[sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
%sC = T.tot_predictions(:,[5 4 7 8 9 10 1]);
sC = T.predictions;
%sC = [T.predictions T.shapeData];
%sC = [T.predictions T.shapeData(:,[2 3 5 6 8 9])];
usC = mean(sC,1);
sC = bsxfun(@minus,sC,usC);
TrainMaster = [];
TestMaster = [];
for u = 1:numel(UQ)
    fidx = find(strcmp(T.genoType,UQ{u}));
    [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
    TrainMaster = [TrainMaster;Train];
    TestMaster = [TestMaster;Test];
end
TrainMaster = find(TrainMaster);
TestMaster = find(TestMaster);
u_sC = mean(sC(TrainMaster,:),1);
Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
u_tC = mean(tC(TrainMaster,:),1);
Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
[mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
Test_mU = Test_sC*mA;
Test_mV = Test_tC*mB;
BETA = mA*inv(mB);
plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;
for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    %UTA = mean(T.tipAngle(fidxw,:),1);
    %STA = std(T.tipAngle(fidxw,:),1,1);
    
    UTA = mean(Y(fidxw,:),1);
    STA = std(Y(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end

% corr with predictions
[RHO PVAL] = corr(mU,sC(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),Y(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end

% regession along vars
for e = 1:size(mU,2)
    CS = linspace(min(mU(:,e)),max(mU(:,e)),10);
    for tm = 1:size(sC,2)
        p = polyfit(mU(:,e),sC(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(tC,2)
        p = polyfit(mV(:,e),tC(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),Y(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
%% look at model factors for only predictions
close all
cv_tipAngle = PCA_BKPROJ(mB',tE,tU);
cv_tipAngle = bsxfun(@minus,cv_tipAngle,tU);
figure;
plot(cv_tipAngle');

%{
t_struct_C = (mV'*tC);
t_struct_C = PCA_BKPROJ(t_struct_C',tE,tU);
t_struct_C = bsxfun(@minus,t_struct_C,tU);
figure;
plot(t_struct_C')
%}

%{
for e = 1:size(mV,2)
    for tm = 1:size(T.tipAngle,2)
        plot(mV(:,e),T.tipAngle(:,tm),'.')
        drawnow
        pause(.1)
    end
end
%}
    
[t_struct_C pval]= corr(mV,T.tipAngle(TrainMaster,:));
figure;
plot(t_struct_C')


s_struct_C = corr(mU,T.predictions(TrainMaster,:));
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

for e = 1:size(mA,2)
    figure;
    bar(mA(:,e))
    L{e} = num2str(e);
    title(['corr_value' num2str(e)]);
end
%% lat parameters alone without hold outs
T = M;
close all
%X = [T.predictions T.shapeData];
X = [T.predictions];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,T.LAT);

[t_struct_C t_pval] = corr(mV,T.LAT);
for e = 1:size(mB,2)
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

[s_struct_C s_pval] = corr(mU,X);
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
end

%{
corr(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3));
figure
plot(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3),'.');
figure;
plot(mU(:,1),mV(:,1),'.');
corr(T.LAT(:,2),T.predictions(:,3));
%}

para = 1;
figure;
kidx = kmeans(mU(:,para),4);
UQ = unique(kidx);
for u = 1:numel(UQ)
    fidx = find(kidx==u);
    uT = mean(T.tipAngle(fidx,:),1);
    uS = std(T.tipAngle(fidx,:),1,1)*numel(fidx).^-.5;
    errorbar(uT,uS);
    hold all
    LEG{u} = num2str(mean(mU(fidx,para)));
end
legend(LEG)
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/lat_core_struct_tip_angle.csv',t_struct_C);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/lat_core_struct_predictions.csv',s_struct_C);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var1.csv',[mU(:,1),mV(:,1)]);
%csvwrite('/mnt/spaldingdata/nate/communications/meetings/talks/lab meeting/13_12_12/scatter_con_var2.csv',[mU(:,2),mV(:,2)]);
%% order NAM parents by cv1 for spec data
UQ = unique(T.genoType);
TROPICAL = [0 0 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 1 0];
%TROPICAL = [0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
TROP = [];
NON_TROP = [];
e_mU = [];
for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    e_mU(u) = mean(mU(fidx));
    %{
    if TROPICAL(u)
        TROP = [TROP;f.tipAngle(fidx,:)];
    else
        NON_TROP = [NON_TROP;f.tipAngle(fidx,:)];
    end
    %}
end
[e_mU sidx] = sort(e_mU);
UQ = UQ(sidx);
close all
bar(e_mU)
%{
e_TROP = mean(TROP,1);
s_TROP = std(TROP,1,1)*size(TROP,1)^-.5;
e_NON_TROP = mean(NON_TROP,1);
s_NON_TROP = std(NON_TROP,1,1)*size(NON_TROP,1)^-.5;
errorbar(e_TROP*180/pi,s_TROP*180/pi,'r');
hold on
errorbar(e_NON_TROP*180/pi,s_NON_TROP*180/pi);
%}
%% predict shape from composition
T = f;
close all
X = [T.predictions];
Y = [T.shapeData];
[mA,mB,mr,mU,mV,mstats] = canoncorr(X,Y);
REL = 3;
[t_struct_C t_pval] = corr(mV,Y);
for e = 1:REL
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
    
end

[s_struct_C s_pval] = corr(mU,X);
for e = 1:REL
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end


for e = 1:size(mA,2)
    figure;
    plot(mU(:,e),mV(:,e),'.')
end
%% master of joy series - hold out 50% genotype of each one
%% 1) random holdout without respect to genotytpe
for tr = 1:100
    T = f;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    %sC = [T.predictions];
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    Test_mU = Test_sC*mA;
    Test_mV = Test_tC*mB;
    UtoV = mU\mV;
    BETA = mA*UtoV*inv(mB);

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    cCore_S_composition_m = corr(mU,T.predictions(TrainMaster,:));
    
    dC(:,:,tr) = differentialCore(Test_sC,Test_tC,mA,mB,cCore_S_composition_m,cCore_T_m,T.predictions(TestMaster,:));
end
%% 1) differntail random holdout with respect to genotytpe
for tr = 1:100
    T = f;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    UQ = unique(f.genoType);
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    
    
    
    %sC = [T.predictions];
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA,mB,mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    Test_mU = Test_sC*mA;
    Test_mV = Test_tC*mB;
    UtoV = mU\mV;
    BETA = mA*UtoV*inv(mB);

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    cCore_S_composition_m = corr(mU,T.predictions(TrainMaster,:));
   
    
    dC(:,:,tr) = differentialCore(Test_sC,Test_tC,mA,mB,cCore_S_composition_m,cCore_T_m,T.predictions(TestMaster,:));
end
%% 2) differntail random holdout with respect to genotytpe
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];
UQ = unique(f.genoType);

for tr = 1:100
    T = M;
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    tC = T.LAT;
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    sC = [T.predictions];
    
    
    %{
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    %}
    %[TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);    
    
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);

    
    if tr == 1
        BASEA = mA;
    else
        SGN(tr) = sign(BASEA(:,e)'*mA(:,e,tr));
        for e = 1:size(BASEA,2)
            if sign(BASEA(:,e)'*mA(:,e,tr)) == -1
                mA(:,e,tr) = -mA(:,e,tr);
                mB(:,e,tr) = -mB(:,e,tr);
                mU(:,e) = -mU(:,e);
                mV(:,e) = -mV(:,e);
            end
        end
        BASEA = mean(mA,3);
    end
    
    
    Test_mU = Test_sC*mA(:,:,tr);
    Test_mV = Test_tC*mB(:,:,tr);;
    UtoV = mU\mV;
    BETA = mA(:,:,tr)*UtoV*inv(mB(:,:,tr));

    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    
    
    
    
    
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    
    
    
    
    
    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    %uK = mean(abs(CORE),3);
    %sK = std(abs(CORE),1,3)
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];
    
    for e = 1:2
        figure(tmp(e));
        bar(uK(e,:))
        hold on
        h = errorbar(uK(e,:),sK(e,:));
        set(h,'LineStyle','none');
        drawnow;
        hold off
    end
end
%% 3) geno typ ehold oout
CORE = [];
tmp(1) = figure;
tmp(2) = figure;
close all
mA = [];
UQ = unique(f.genoType);
PRE = [];
ACC = [];
T = f;
for tr = 1:numel(UQ);
    
    %initT = T.tipAngle(:,1);
    %T.tipAngle = bsxfun(@minus,T.tipAngle,initT);
    perDraw = .5;
    [tS tC tU tE tL tERR tLAM] = PCA_FIT_FULL(T.tipAngle,3);
    tC = T.LAT;
    [sS sC sU sE sL sERR sLAM] = PCA_FIT_FULL(T.specData,50);
    sC = [T.predictions];
    
    
    %{
    TrainMaster = [];
    TestMaster = [];
    for u = 1:numel(UQ)
        fidx = find(strcmp(T.genoType,UQ{u}));
        [Train, Test] = crossvalind('HoldOut', size(fidx,1),perDraw);
        TrainMaster = [TrainMaster;Train];
        TestMaster = [TestMaster;Test];
    end
    %}
    %[TrainMaster, TestMaster] = crossvalind('HoldOut', size(T.tipAngle,1),perDraw);    
    TrainMaster = ~strcmp(T.genoType,UQ{tr});
    TestMaster = strcmp(T.genoType,UQ{tr});
    TrainMaster = find(TrainMaster);
    TestMaster = find(TestMaster);
    u_sC = mean(sC(TrainMaster,:),1);
    Train_sC = bsxfun(@minus,sC(TrainMaster,:),u_sC);
    Test_sC = bsxfun(@minus,sC(TestMaster,:),u_sC);
    u_tC = mean(tC(TrainMaster,:),1);
    Train_tC = bsxfun(@minus,tC(TrainMaster,:),u_tC);
    Test_tC = bsxfun(@minus,tC(TestMaster,:),u_tC);
    [mA(:,:,tr),mB(:,:,tr),mr,mU,mV,mstats] = canoncorr(Train_sC,Train_tC);
    
    
    [gmA,gmB,gmr,gmU,gmV,gmstats] = canoncorr(sC,tC);
    gCORE = corr(gmU,T.predictions);
    
    if tr == 1
        BASEA = mA;
    else
        SGN(tr) = sign(BASEA(:,e)'*mA(:,e,tr));
        for e = 1:size(BASEA,2)
            if sign(BASEA(:,e)'*mA(:,e,tr)) == -1
                mA(:,e,tr) = -mA(:,e,tr);
                mB(:,e,tr) = -mB(:,e,tr);
                mU(:,e) = -mU(:,e);
                mV(:,e) = -mV(:,e);
            end
        end
        BASEA = mean(mA,3);
    end
    
    
    Test_mU = Train_sC*mA(:,:,tr);
    Test_mV = Test_tC*mB(:,:,tr);;
    UtoV = mU\mV;
    BETA = mA(:,:,tr)*UtoV*inv(mB(:,:,tr));

    
    
    PRE = [PRE;Test_sC*mA(:,:,tr)];
    ACC = [ACC;Test_tC*mB(:,:,tr)]
    
    
    
    
    
    
    
    Core_T_m = mV\tC(TrainMaster,:);
    Core_S_m = mU\tS(TrainMaster,:);
    
    
    
    
    
    Core_S_composition_m = mU\T.predictions(TrainMaster,:);

    
    
    
    
    
    cCore_T_m = corr(mV,tC(TrainMaster,:));
    cCore_T_m = corr(mV,T.tipAngle(TrainMaster,:));
    cCore_S_m = corr(mU,tS(TrainMaster,:));
    CORE(:,:,tr) = corr(mU,T.predictions(TrainMaster,:));
    
    uCORE = squeeze(mean(CORE,3));
    for e = 1:size(CORE,1)
        if sign(uCORE(e,:)*squeeze(CORE(e,:,tr))') == -1
            CORE(e,:,tr) = -CORE(e,:,tr);
        end
    end
    
    %uK = mean(abs(CORE),3);
    %sK = std(abs(CORE),1,3)
    uK = mean((CORE),3);
    sK = std((CORE),1,3);
    K = [sK.*uK.^-1];
    
    for e = 1:2
        figure(tmp(e));
        bar(uK(e,:))
        hold on
        h = errorbar(uK(e,:),sK(e,:));
        set(h,'LineStyle','none');
        drawnow;
        hold off
    end
end

%%
for e = 1:size(CORE,3)
    for cc = 1:size(CORE,1)
        if sign(squeeze(CORE(cc,:,e))*gCORE(cc,:)') == -1
            CORE(cc,:,e) = -CORE(cc,:,e);
        end
    end
end
K = abs(bsxfun(@minus,CORE,gCORE));
K = squeeze(K);
K = squeeze(sum(K,2));
sCORE = sum(gCORE,2);
K = bsxfun(@times,K,sCORE.^-1);
[s sidx] = sort(K(1,:));
sUQ = UQ(sidx);
figure
bar(s)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))
sortCORE = CORE(:,:,sidx);




for u = 1:numel(UQ)
    fidx = strcmp(T.genoType,UQ{u});
    mv(u,:) = sum(gmU(fidx,:));
end
[s sidx] = sort(mv(:,1));
sUQ = UQ(sidx);
figure
bar(s)
set(gca,'XTickLabel',sUQ,'XTick',1:numel(sUQ))

figure;



%% JUNK








for e = 1:size(tC,2)
    CS = linspace(min(tC(TrainMaster,e)),max(tC(TrainMaster,e)),10);
    SIM = Core_T(e,:)'*CS;
    predict_tipAngle = PCA_BKPROJ(SIM',tE,tU);
    figure;plot(predict_tipAngle')
end


plot(mU(:,1),mV(:,1),'.');
hold on
plot(Test_mU(:,1),Test_mV(:,1),'r.');
figure;


Prediction_tC = Test_sC*mA;
delta = Prediction_tC - Test_tC;
delta = sum(delta.*delta,2).^.5;

for e = 1:size(Prediction_tC,2)
    figure;
    plot(Prediction_tC(:,e),Test_tC(:,e),'.')
end

plot(Prediction_tC(:,1),Test_tC(:,1),'.')








for u = 1:numel(UQ)
    fidxw = find(strcmp(T.genoType,UQ{u}));
    fidx = intersect(fidxw,TestMaster);
    predict = bsxfun(@minus,sC(fidx,:),u_sC)*BETA;
    predict_tipAngle = PCA_BKPROJ(predict,tE,tU);
    pUTA = mean(predict_tipAngle,1);
    pSTA = std(predict_tipAngle,1,1);
    
    UTA = mean(T.tipAngle(fidxw,:),1);
    STA = std(T.tipAngle(fidxw,:),1,1);
    
    errorbar(pUTA,pSTA,'r');
    hold on
    errorbar(UTA,STA,'k');
    hold off
    pause(.1);
end


% corr with predictions
[RHO PVAL] = corr(mU,T.predictions(TrainMaster,:));

% regession along vars
for e = 1:size(mV,2)
    CS = linspace(min(mV(:,e)),max(mV(:,e)),10);
    for tm = 1:size(T.tipAngle,2)
        p = polyfit(mV(:,e),T.tipAngle(TrainMaster,tm),1);
        RVAL(tm,:,e) = corr(mV(:,e),T.tipAngle(TrainMaster,tm));
        SIMU(tm,:,e) = polyval(p,CS);
    end
end
 



%%
T = f;
close all

[mA,mB,mr,mU,mV,mstats] = canoncorr(T.predictions,T.shapeData);

[t_struct_C t_pval] = corr(mV,T.LAT);
for e = 1:size(mB,2)
    figure;
    bar(t_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end

[s_struct_C s_pval] = corr(mU,T.predictions);
for e = 1:size(mA,2)
    figure;
    bar(s_struct_C(e,:))
    title(['struct_core' num2str(e)]);
end
corr(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3));
figure
plot(mean(T.tipAngle(:,end-10:end),2),T.predictions(:,3),'.');
figure;
plot(mU(:,1),mV(:,1),'.');
