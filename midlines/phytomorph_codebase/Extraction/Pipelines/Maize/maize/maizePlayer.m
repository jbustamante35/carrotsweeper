%% maize player
close all
sPath = '/mnt/scratch5/maizeContours/';
%sPath = '/mnt/scratch5/maizeContours_withmidlines/';
%%% scan for new images
FileList = {};
FileExt = {'mat'};
verbose = 1;
SET = gdig(sPath,FileList,FileExt,verbose);
L = 1;
%%
for e = 10%:numel(SET)
    try
        %figure;
        root = load(SET{e});
        root = root.obj;
        for loop = 1:L
            root.plot();
            %angle = root.getTipAngle();
            %length = root.getMidlineLength();
        end
        %{
        %figure;
        root = load(strrep(SET{e},'maizeContours2','maizeContours'));
        root = root.obj;
        for loop = 1:L
            root.plot();
        end
        %}
    catch ME
        
    end
end
%% extract angles alone
close all
sPath = '/mnt/scratch5/maizeContours2/';
sPath = '/mnt/scratch5/maizeContours_withmidlines/';
%%% scan for new images
FileList = {};
FileExt = {'mat'};
verbose = 1;
SET = gdig(sPath,FileList,FileExt,verbose);
L = 1;
for e = 1:numel(SET)
    try
        %figure;
        root = load(SET{e});
        root = root.obj;
        angle(e,:) = root.getTipAngle();
        length(e,:) = root.getMidlineLength();
    catch ME
        
    end
end
%% 
dA = abs(diff(angle,1,1))