%%
user = 'nhaase';
resultsPath = ['/mnt/spaldingdata/nate/mirror_images/maizeData/' user '/return/earData/output/'];
% dig for images
FilePath = strrep(inPath,'#USER#',user);    
FileList = {};
FileExt = {'mat'};
verbose = 1;
[FileList] = gdig(resultsPath,FileList,FileExt,verbose);
%%
N = [];
O = [];
for e = 1:numel(FileList)
        try
        RAD = 1200:25:1600;
        data = load(FileList{e});
        data.KernelLength;
        newT = [];
        for ear = 1:3
            for r = 1:numel(RAD)
                dR = RAD(r);
                ret = data.FT.G{ear}{r};
                ufT = ret(1:dR);
                h = fspecial('average',[5 1]);
                ufT = imfilter(ufT,h);
                [newT(ear,r) f] = findT(ufT,2*dR+1);
            end
        end    
        oldT = mean(data.KernelLength,2);
        newT = mean(newT,2);
        if any(newT > 700)
            hello = 1
        end
        N = [N;newT(:)];
        O = [O;oldT(:)];
        e
        catch
            
        end
end