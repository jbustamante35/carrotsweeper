function [] = singlePHstack(imageStack,oPath)
    mkdir(oPath);

    norLEN = 6;
    norScale = linspace(0,norLEN,norLEN);
    sampleAlong = 450;
    for e = 1:numel(imageStack)
        tic
        fileName = imageStack{e};
        tdb{e} = extractBoundary(fileName,11);        
        %tdb{e} = smoothCurve(tdb{e});
        [tC{e}] = getTip(tdb{e},31);
        [tN{e}] = generateNormalField(tdb{e},7);
        [tS{e} tM{e} tT{e} tB{e}] = sampleBoundary(fileName,tdb{e},tN{e},norScale,tC{e},sampleAlong,0);
        toc
    end


    S = [];
    for e = 1:numel(imageStack)
        S = cat(3,S,tS{e});
    end

    S = permute(S,[1 3 2]);    
    ratio = S(:,:,2).*S(:,:,1).^-1;
    R = ratio;
    h = fspecial('average',[5 5]);
    csvwrite([oPath '_ratio.csv'],R);
    %mesh(imfilter(R,h,'replicate'));
    %view([0 90]);
    %drawnow


    for e = 1:numel(imageStack)

        I = imread(imageStack{e});
        I = I(:,:,1:3);
        nSKIP = 10;
        image(I);
        hold on        
        ED = [];
        for p = 1:size(tdb{e}{1},1)
            str = tdb{e}{1}(p,:);
            ed = str + norScale(end)*tN{e}(p,:);
            if mod(p,nSKIP) == 0
                plot([str(1) ed(1)],[str(2) ed(2)],'r')
            end
            ED = [ED;ed];
        end



        plot(ED(:,1),ED(:,2),'r');
        plot(tdb{e}{1}(:,1),tdb{e}{1}(:,2),'r');

        plot(ED(tT{e}(1):tC{e},1),ED(tT{e}(1):tC{e},2),'c');
        plot(tdb{e}{1}(tT{e}(1):tC{e},1),tdb{e}{1}(tT{e}(1):tC{e},2),'c');


        plot(ED(tC{e}:tB{e}(end),1),ED(tC{e}:tB{e}(end),2),'b');
        plot(tdb{e}{1}(tC{e}:tB{e}(end),1),tdb{e}{1}(tC{e}:tB{e}(end),2),'b');



        plot(tdb{e}{1}(tC{e},1),tdb{e}{1}(tC{e},2),'g*');

        drawnow
        hold off
        saveas(gca,[oPath 'frame_' num2str(e) '.tif']);
        close all
    end

    here =1;
end

%{

    FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/pH-measurements/';
    FilePath = '/mnt/scratch3/users/monshausenlab/kinematics/Gaby/test4/';
    FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/';
    FilePath = '/home/nate/Downloads/Fusicoccin/';
    FilePath = '/home/nate/Downloads/newGabyPH/Fusi/';
    FilePath = '//home/nate/Downloads/newGabyPH/Gravi/';
    %FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Gravitropism/Gaby/';
    %FilePath = '/mnt/snapper/nate/mirror_images/Surface pH/Auxin treatment/Gaby/';
    FileList = {};
    FileExt = {'tif','TIF'};
    verbose = 1;
    FileList = sdig(FilePath,FileList,FileExt,verbose);

    singlePHstack(FileList{1})

%}