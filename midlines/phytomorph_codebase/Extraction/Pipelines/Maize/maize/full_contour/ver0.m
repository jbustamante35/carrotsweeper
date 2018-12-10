FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
FilePath = '/mnt/spaldingdata/Takeshi/B73 vs NC350 (20130605-)/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% read in the first image from each stack
% preallocate
I = imread(SET{1}{1});
S = zeros([size(I) numel(SET)]);

for e = 1:numel(SET)
    % 
    %S(:,:,e) = imread(SET{e}{1});
    %S(:,:,e) = S(:,:,e)/255;
    
    I = double(imread(SET{e}{1}));
    I = I/255;
    
    
    % get level contours
    curves = getLevelContours(I,10);
    % remove non closed curves
    curves = selectClosedCurves(curves);
    % remove small curves
    ridx = [curves.length] < 200;
    curves(ridx) = [];
    % sample curve bank
    curves = sampleCurveBank(I,curves);
    
    D{e}.curves = curves;
    
    e
    numel(SET)
end
%%
for e = 1:numel(SET)
    I = double(imread(SET{e}{1}))/255;
    T = 30;
    R = 50;
    [THETA RAD] = ndgrid(linspace(-pi,pi,T),linspace(0,30,R));
    NH = [RAD(:).*cos(THETA(:)) RAD(:).*sin(THETA(:))]';
    sam{e} = bugEye(I,NH);
    e
    numel(SET)
end