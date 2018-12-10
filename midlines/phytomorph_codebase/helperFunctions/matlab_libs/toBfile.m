function [tList] = toBfile(fileList)
    
    import java.util.ArrayList;    
    import phytoG.locked.Bpersist.Bfs.abstractions.*;
    tList = ArrayList();
    if isa(fileList,'cell')
        for f = 1:numel(fileList)
            tmp = Buri_file(fileList{f});
            tList.add(tmp);
        end
    elseif isa(fileList,'imageStack')
        for f = 1:numel(fileList)
            tmp = Buri_file(fileList{f}.fileName);
            tList.add(tmp);
        end
    end
        
end