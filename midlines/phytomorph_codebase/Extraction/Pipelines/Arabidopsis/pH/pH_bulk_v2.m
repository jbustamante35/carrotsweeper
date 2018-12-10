%% 1) look for images
FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/pH-measurements/';
FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/test4/';
FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/';
FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Gravitropism/Gaby/';
FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Auxin treatment/Gaby/';
FileList = {};
FileExt = {'tif','TIF'};
verbose = 1;
FileList = sdig(FilePath,FileList,FileExt,verbose);
%% 2) watch movies
SKIP = 10;
for s = 1:numel(FileList)
    for e = 1:SKIP:numel(FileList{s})
        I = imread(FileList{s}{e});
        imshow(I,[]);
        drawnow
    end
end
%% 2.5) remove 1122 data
moniker = '1122';
rmidx = [];
for e = 1:numel(FileList)
    [p,n,ext] = fileparts(FileList{e}{1});
    if ~isempty(strfind(p,moniker))
       rmidx = [rmidx e];
    end
end
FileList(rmidx) = [];

%% find cngc
moniker = 'cngc';
subSet = [];
for e = 1:numel(FileList)
    if ~isempty(strfind(FileList{e}{1},moniker))
       subSet = [subSet e];
    end
end
%% 2.6) calibrate
close all
phD = readtext('/home/nate/Downloads/pH calibration - 2014-12 (center).csv');
PH = cell2mat(phD(:,1));
uPH = mean(PH);
RATIO = cell2mat(phD(:,2));
uR = mean(RATIO);
PH = PH - uPH
RATIO = RATIO - uR;
sigfunc = @(A, x)A(5) + (A(1) ./ (A(2) + exp(-A(3)*(x-A(4)))));
A0 = [1 1 1 1 0]; %// Initial values fed into the iterative algorithm
fit = nlinfit(PH, RATIO, sigfunc, A0);
iPH = linspace(4-uPH,8-uPH,10000);
Y = sigfunc(fit,iPH);
plot(PH+uPH,RATIO+uR,'b*');
hold on
plot(iPH+uPH,Y+uR,'r');
iPH = iPH + uPH;
Y = Y + uR;

%% 3) extract root
disp = 0;
for s = 1:numel(FileList)
    s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear samp    
    norScale = linspace(-10,30,40);
    sampleAlong = 400;
    parfor e = 1:numel(FileList{s})        
        fileName = FileList{s}{e};
        dB{e} = extractBoundary(fileName,11);        
        dB{e} = smoothCurve(dB{e});
        [tipC{e}] = getTip(dB{e},31);
        [NOR{e}] = generateNormalField(dB{e},7);
        [samp{e} mask{e} topIdx{e} botIdx{e}] = sampleBoundary(fileName,dB{e},NOR{e},norScale,tipC{e},sampleAlong,0);        
        numel(FileList{s})
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display heat map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = [];
    for e = 1:numel(FileList{s})
        S = cat(3,S,samp{e});
    end
    S = permute(S,[1 3 2]);    
    ratio = S(:,:,2).*S(:,:,1).^-1;
    phV = [];
    for e = 1:numel(ratio)
        delta = abs(ratio(e) - Y);
        [~,midx] = min(delta);
        phV(e) = iPH(midx);
    end
    phV = reshape(phV,size(ratio));
    
    
    PHDATA{s} = phV;
    
    if disp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display movie with contour
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nSKIP = 10;
        for e = 1:numel(FileList{s})
            I = imread(FileList{s}{e});
            imshow(I,[]);
            hold on        
            ED = [];
            for p = 1:size(dB{e}{1},1)
                str = dB{e}{1}(p,:);
                ed = str + norScale(end)*NOR{e}(p,:);
                if mod(p,nSKIP) == 0
                    plot([str(2) ed(2)],[str(1) ed(1)],'r')
                end
                ED = [ED;ed];
            end



            plot(ED(:,2),ED(:,1),'r');
            plot(dB{e}{1}(:,2),dB{e}{1}(:,1),'r');

            plot(ED(topIdx{e}(1):tipC{e},2),ED(topIdx{e}(1):tipC{e},1),'c');
            plot(dB{e}{1}(topIdx{e}(1):tipC{e},2),dB{e}{1}(topIdx{e}(1):tipC{e},1),'c');


            plot(ED(tipC{e}:botIdx{e}(end),2),ED(tipC{e}:botIdx{e}(end),1),'b');
            plot(dB{e}{1}(tipC{e}:botIdx{e}(end),2),dB{e}{1}(tipC{e}:botIdx{e}(end),1),'b');



            plot(dB{e}{1}(tipC{e},2),dB{e}{1}(tipC{e},1),'g*');

            drawnow
            hold off
        end
    end
        
end
%%
SMOOTH = 5;
tipK = 151;
bulkC = 7;
sampleLength = 400;

for e = 1:numel(FileList)
    try
        disp = 0;
        rootCurve = {};
        parfor i = 1:numel(FileList{e})
            I = imread(FileList{e}{i});
            I = flipdim(I,1);
            rootCurve{i} = isolateRoot2(I,bulkC);
            if disp
                imshow(I,[]);
                hold on;
                plot(rootCurve{i}(1,:),rootCurve{i}(2,:),'r');
                drawnow
            end
             i
            numel(FileList{e})
        end
        [tipC] = getTip(rootCurve,tipK);
        [SAM(e).data] = sampleRGB(FileList{e},rootCurve,tipC,sampleLength,SMOOTH);
        SAM(e).data{1} = flipdim(SAM(e).data{1},1);
        [dP] = measureDisplacement(rootCurve,tipC);
        
        hold on
        for i = 1:numel(FileList{e})
           I = imread(FileList{e}{i});
           I = flipdim(I,1);
           imshow(I,[]);
           hold on
           plot(rootCurve{i}(:,2),rootCurve{i}(:,1));
           plot(rootCurve{i}((tipC(i)-sampleLength):tipC(i),2),rootCurve{i}((tipC(i)-sampleLength):tipC(i),1),'k','LineWidth',3);
           plot(rootCurve{i}((tipC(i):tipC(i)+sampleLength),2),rootCurve{i}(tipC(i):(tipC(i)+sampleLength),1),'m','LineWidth',3);
           hold on
           plot(rootCurve{i}(tipC(i),2),rootCurve{i}(tipC(i),1),'r*');
           hold off
           title(i)
           drawnow
        end
    
        figure;
        mesh(SAM(e).data{1}(:,:,3));
        view([0 90]);
        figure;
        mesh(SAM(e).data{2}(:,:,3));
        view([0 90]);
    
        
        waitforbuttonpress
    catch ME
        ME
    end
   
    
    
    
end
%%