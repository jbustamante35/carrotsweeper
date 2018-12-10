function [] = addMetaToFileList(FileList,algo,ver)
    mainline = 'imeta add -d #dataObject# #key# #value# #nvalue#';
    for e = 1:numel(FileList)
        CMD = strrep(mainline,'#dataObject#',FileList{e});
        CMD = strrep(CMD,'#key#','algo');
        CMD = strrep(CMD,'#value#',[algo '.' ver]);
        
    end
end