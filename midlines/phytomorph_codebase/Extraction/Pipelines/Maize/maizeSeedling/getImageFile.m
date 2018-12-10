function [imageFile] = getImageFile(FileList,PlotNumber,PictureDay)
    key1 = 'Plot';
    key2 = 'PictureDay';
    imageFile = '';
    for e = 1:numel(FileList)
        v1 = keyLookup(FileList{e},key1);
        if strcmp(v1,PlotNumber)
            v2 = keyLookup(FileList{e},key2);
            if strcmp(v2,PictureDay)
                [p,nm,ext] = fileparts(FileList{e});
                if strcmp(ext,'.txt')
                    CMD = ['cat "' FileList{e} '"'];
                    [status,imageFile] = system(CMD);
                    fidx = strfind(imageFile,'#');
                    if ~isempty(fidx)
                        imageFile(fidx(1):end) = [];
                    end
                    imageFile = strrep(imageFile,char(10),'');
                    return
                end
            end
        end
    end
end