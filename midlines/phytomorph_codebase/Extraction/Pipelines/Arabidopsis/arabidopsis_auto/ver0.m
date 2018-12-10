%% find files and generate rand number
pth = '/mnt/scratch3/users/nmiller/phytoMorph/outPort/'
dataFiles = gdig(pth,{},{'mat'},1);
numel(dataFiles);
%% domain disk and square
Tr = tracer();
Tr.setXsz([81 81]);
Tr.setYsz([41 2]);
Tr.setinitComp(20);
Tr.setGroupN(9);
Tr.setmodelCompX(8);
Tr.setmodelCompY(5);

p.type = 'disk';
p.value{1} = [0 30 200];
p.value{2} = [-pi pi 100];
tD = phytoAdomain(p);
tD.generateDomain();

p.type = 'box';
p.value{1} = [-40 40 81];
p.value{2} = [-40 40 81];
tD = phytoAdomain(p);
tD.generateDomain();
%%
for e = 1:10
    load(dataFiles{e})

    ridx = randi(300,2,1);
    for i = 1:numel(ridx)
        
        close all;
        rootN = ridx(i);
        
        I = imread(d{rootN}{1}.image.fileName);
        [I direc] = handleFLIP(I,d{rootN}{1}.flipD);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        try
            % get tangent space
            TM = d{rootN}{1}.midline_1.generateTangentSpace(15);
            patchS = TM.sample(I,tD);

            % lift curve
            gamma = liftCurve(d{rootN}{1}.midline_1,40,TM);
            gamma = reshape(gamma,[size(gamma,1) size(gamma,2)*size(gamma,3)]);
        catch ME
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add to X and Y
        Tr.addXY(patchS.d,gamma(1:size(patchS.d,1),:,:));
        
        % show 
        h = figure;
        imshow(I);
        hold on
        d{rootN}{1}.midline_1.view(h,[],[]);
        TM.view(h,[],[])
        axis equal
        
        drawnow
    end
end
%%
Tr.learn();
%%
Tr.viewManifold();
%% patch

h = figure;
patchS.view(h,[],[]);

%% look at curves

close all

Tr.viewCoandDo();
close all




