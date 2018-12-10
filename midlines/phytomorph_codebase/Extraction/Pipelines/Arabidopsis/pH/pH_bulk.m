FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/pH-measurements/';
FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/test4/';
FileList = {};
FileExt = {'tif','TIF'};
verbose = 1;
FileList = sdig(FilePath,FileList,FileExt,verbose);

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