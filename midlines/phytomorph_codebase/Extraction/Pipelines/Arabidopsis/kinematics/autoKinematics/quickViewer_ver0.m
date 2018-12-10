%% load from condor
FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/monshausenlab/return/straightKinematicsWP/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/arabidopsisData/trieupham10/return/straightKinematicsQ1/';
FileList = {};
FileExt = {'mat'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
%% 
oPath = '/home/nate/Downloads/forGaby2/';
mkdir oPath
close all
h1 = figure;
h2 = figure;
h3 = figure;
G = {};
LEG = {};
for e = 1:numel(FileList)
    close all
    
    F = load(FileList{e},'F','nm','pointList');
    
    if 1
       
    %if ~isempty(strfind(F.nm,'Nathan'))%~isempty(strfind(F.nm,'highres')) %& ~isempty(strfind(F.nm,'WT')) %& ~isempty(strfind(F.nm,'20160527'))  
    %if ~isempty(strfind(F.nm,'Hypo and Hyper'))
    if ~isempty(strfind(F.nm,'16')) & ~isempty(strfind(lower(F.nm),'apy1oe'))
    F.nm
    ridx = strfind(F.nm,'#');
    F.nm(ridx(1):end) = [];
    I = imread(F.nm);
        
        [p,n,ext] = fileparts(F.nm);
        fidx = strfind(p,filesep);
        n = p(fidx(end)+1:end);
        F.nm
        %{
        FilePath = [p filesep];
        FileList = {};
        FileExt = {'TIF','tif','png'};
        FileList = gdig(FilePath,FileList,FileExt,1);
        
        %}
        
        %{
        [J(:,:,1) BOX] = imcrop(I);
        for e = 1:numel(FileList)
            I = imread(FileList{e});
            J(:,:,e) = imcrop(I,BOX);
        end
        
        
        
        
        
        for loop = 1:1
            for e = 1:size(J,3)
                imshow(J(:,:,e),[])
                drawnow
            end
        end    
        %}
        %{
        figure(h1);
        imshow(I,[]);
        hold on
        plot(F.pointList(:,2,1),F.pointList(:,1,1),'r')
        for e = 1:10
            imshow(I,[])
            hold on
            plot(F.pointList(:,2,e),F.pointList(:,1,e),'r.')
            drawnow
            waitforbuttonpress
            hold off
        end
        
        waitforbuttonpress
        %}
        imshow(I,[]);
        hold on
        plot(F.pointList(1:end,2,2),F.pointList(1:end,1,2),'r','LineWidth',3)
        %[domain,vP] = measureVelocityProfile(F.pointList(20:end,:,:),1,1,10);
        %F.F = vP;
        sr = extractStrain(F.F(:,:,end),[21 3]);
        figure(h2);
        K = interp2(sr(:,:,2),3);
        mesh(K);
        view([0 90]);
        axis([1 size(K,2) 1 size(K,1)]);
        colorbar
        %saveas(gca,[oPath n '-colorbar.tif']);
        
        figure(h3);
        plot(mean(sr(:,:,2),2))
        %saveas(gca,[oPath n '-colorbar.tif']);
        
        [p,n,ex] = fileparts(F.nm);
        fidx = strfind(p,filesep);
        T = p(fidx(end)+1:end);
        T = strrep(T,'_','--')
        title(T);
        waitforbuttonpress
        %{
        %waitforbuttonpress
        button = questdlg('Keep?');
        if strcmp(button,'Yes')
            G{end+1} = mean(sr(:,:,2),2);
            LEG{end+1} = T;
        end
        %}
    end
    end
end
%% 
close all
figure;
LEGG = {};
for e = 1:numel(G)
    try
        K(:,e) = G{e}(1:80);
        plot(G{e}(1:80));
        LEGG{e} = LEG{e};
        hold all
    catch
    end
end
legend(LEGG)
%% 
close all
searchFor = {'wswt','apy1','apy2'};
CL = [];
for e = 1:numel(LEG)
    for k = 1:numel(searchFor)
        if ~isempty(strfind(lower(LEGG{e}),searchFor{k}))
            CL(e) = k;
        end
    end
end
UQ = unique(CL);
for u = 1:numel(UQ)
    fidx = find(CL==UQ(u));
    U = mean(K(:,fidx),2);
    S = std(K(:,fidx),1,2)*numel(fidx)^-.5;
    errorbar(U,S);
    hold all
end
legend(searchFor)
