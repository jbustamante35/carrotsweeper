colGroups = {'garf0012' 'nhaase' 'gxe' 'ruairidh'};
basePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/$user/return/';
for e = 1:numel(colGroups)
   PTH{e} = strrep(basePath,'$user',colGroups{e})
end
%%
CMD = 'find $1';
for e = 1:numel(PTH)
    cmd = strrep(CMD,'$1',PTH{e});
    [status,results{e}] = system(cmd);
end
%%
%% parse out the file types
for r = 1:numel(results)
    nidx = strfind(results{r},char(10));
    str = 1;
    words = {};
    for e = 1:numel(nidx)
        stp = nidx(e);
        fileName = results(str:stp-1);
        [pth nm ext] = fileparts(fileName);
        if strcmp(ext,'.mat')
            fl{r}{end+1} = fileName;
            str = stp + 1;
        end
    end
end