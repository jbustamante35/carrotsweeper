FilePath = '/mnt/snapper/nate/sorghumDrone/images/';
subFolder = 'fourthFlight_energy_am_1st/';
FilePath = [FilePath subFolder];
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% look at 
%% view "movie"
for e = 1:10:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%%
close all
N = 400;
I1 = imread(FileList{N});
I2 = imread(FileList{N+1});
for e = 1:100
    if mod(e,2) == 0
        imshow(I1,[]);
        drawnow
    else
        imshow(I2,[]);
        drawnow
    end
end
%%
rPath = ['/mnt/snapper/nate/sorghumDrone/cimages/' subFolder];
mkdir(rPath);
parfor e = 1:numel(FileList)
    tic
    I = double(imread(FileList{e}))*(2^16-1)^-1;
    [pth,n,ext] = fileparts(FileList{e});
    U = mean(I(:));
    Io = I;
    m = min(I(:));
    I = I - m;
    M = max(I(:));
    I = I / M;
    %I = I / M;
    pd = fitdist(I(:),'rayleigh');
    for w = 8:16
        tmpP = [rPath 'windowLevel_' num2str(w) '/'];
        mkdir(tmpP)
        tmpF= [tmpP n '.tif'];
        J = adapthisteq(I,'NumTiles',[w w],'NBins',2^16);
        J = J * M;
        J = J + m;
        J = J * 2^16;
        J = uint16(J);
        imwrite(J,tmpF);
        %J1 = adapthisteq(I,'NumTiles',[w w],'Distribution','rayleigh','Alpha',pd.B,'NBins',2^16);
        %J1 = J1 * M;
        %J1 = J1 + m;
    end
    e
    toc
end