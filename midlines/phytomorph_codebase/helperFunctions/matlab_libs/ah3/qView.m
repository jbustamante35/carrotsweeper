function [] = qView(SAM,F)
    F.D = reshape(SAM.dVec,[SAM.sz size(SAM.dVec,2)]);    
    F.op = 1;
    F = stateChange(F);
    data = F.D;
    newSZ = size(data);
    data = reshape(data,[prod(newSZ(1:end-1)) newSZ(end)]);
    UQ = unique(SAM.gVec(:,1));
    for loop = 1:50
        for gp = 1:numel(UQ)
            gidx = (SAM.gVec(:,1) == gp);
            sd = data(:,gidx);
            sp = gprobData(sd',F.model);
            imgFile = SAM.fn{gp};
            I = imread(imgFile);
            pt = SAM.pt(gidx,:);
            e = maizeKN(char(imgFile.getFullFileName()));
            [J sidx] = sort(sp,'descend');
            subidx = sidx(1:e);
            imshow(I,[]);
            hold on
            plot(pt(:,1),pt(:,2),'go');
            plot(pt(subidx,1),pt(subidx,2),'ro','MarkerSize',10);
            pause(1)
            hold off
        end
    end
end